# T7/Illumina 测序数据自动化处理流程

本说明文档给出了一套适用于华大 **T7** 平台（FASTQ 格式与 Illumina 兼容）的变异检测流程。流程覆盖了从原始下机数据（`*.fastq.gz`）到经过质量控制、比对、变异检测与过滤的 `VCF` 文件的完整步骤，同时提供了可执行脚本示例与常见可选操作。

## 1. 数据准备

1. 将 `R1/R2` 原始测序文件（`sample_R1.fastq.gz`, `sample_R2.fastq.gz`）放入同一目录。
2. 准备参考基因组 `ref-202509.fna` 并确保其为 `FASTA` 格式。若文件来自 NCBI 或华大自建数据库，需要与变异数据库（如 dbSNP）保持版本一致。
3. 可选：准备用于 BQSR 的已知变异位点（`known_sites.vcf.gz`，可包含 dbSNP、1000G 资源）。
4. 确认软件已安装并加入 `PATH`：详见项目根目录的 [`requirements.txt`](../requirements.txt)。

> **T7 数据兼容说明**：T7 平台输出的 FASTQ 格式与 Illumina 兼容，可直接使用 `bwa`, `fastp` 等主流工具处理。唯一需要关注的是 read group 信息，可结合测序批次自行设置。

## 2. 快速开始：使用脚本一键运行

`scripts/run_pipeline.sh` 封装了各步骤，建议在准备好环境后直接调用。

```bash
bash scripts/run_pipeline.sh \
  -1 /path/to/sample_R1.fastq.gz \
  -2 /path/to/sample_R2.fastq.gz \
  -r /path/to/ref-202509.fna \
  -o results/sample01 \
  -s sample01 \
  -k /path/to/known_sites.vcf.gz \
  -t 32 \
  -p /home/wslll/13_anti_disease_natural/0_biosoftware/picard.jar \
  -g '@RG\tID:lane1\tSM:sample01\tPL:BGISEQ'
```

关键参数说明：

- `-t`：多线程并行数，根据服务器核心数调整。
- `-p`：Picard jar 路径，若已在 `PATH` 中可省略。
- `-g`：可选的 read group 描述，建议在 T7 数据中将 `PL` 设置为 `BGISEQ` 以便下游分析区分平台。

脚本会自动执行以下步骤：

1. `FastQC` 原始 reads。
2. `fastp` 过滤、剪切并生成报告，再次 `FastQC` 检查。
3. 可选 `MultiQC` 汇总质控报告。
4. 构建/检查参考基因组索引（`bwa index`、`samtools faidx`、`gatk CreateSequenceDictionary`）。
5. `bwa mem` 比对 + `samtools sort` 排序，输出 `sorted.bam`。
6. `Picard MarkDuplicates` 标记重复并生成指标。
7. 可选 `GATK BaseRecalibrator/ApplyBQSR` 进行碱基质量再校正。
8. `GATK HaplotypeCaller` 生成 `gVCF`；`GATK GenotypeGVCFs` 生成原始 `VCF`。
9. `GATK VariantFiltration` + `bcftools view -f PASS` 获得过滤后的 `VCF`。
10. `vcftools` 生成深度/质量统计表。

所有中间结果与日志均在 `-o` 指定的输出目录下分类保存：

- `qc/`：FastQC 报告
- `trimmed/`：过滤后的 fastq 及 fastp 报告
- `alignment/`：比对 BAM、重复标记结果、Flagstat 指标
- `gatk/`：BQSR 中间文件
- `variants/`：gVCF、原始/过滤后 VCF 以及统计结果
- `report/`：MultiQC 汇总报告（若安装）

## 3. 手动执行（可替代脚本）

若需逐步调试，可以手动按步骤运行：

### 3.1 质控与清洗
```bash
fastqc -t 16 -o qc/raw sample_R1.fastq.gz sample_R2.fastq.gz
fastp \
  --thread 16 \
  --detect_adapter_for_pe \
  --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 50 \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o trimmed/sample_R1.trimmed.fastq.gz \
  -O trimmed/sample_R2.trimmed.fastq.gz \
  --html trimmed/sample.fastp.html --json trimmed/sample.fastp.json
fastqc -t 16 -o qc/trimmed trimmed/sample_R1.trimmed.fastq.gz trimmed/sample_R2.trimmed.fastq.gz
```

### 3.2 参考基因组索引
```bash
bwa index ref-202509.fna
samtools faidx ref-202509.fna
gatk CreateSequenceDictionary -R ref-202509.fna -O ref-202509.dict
```

### 3.3 比对与 BAM 处理
```bash
bwa mem -t 32 ref-202509.fna trimmed/sample_R1.trimmed.fastq.gz trimmed/sample_R2.trimmed.fastq.gz \
  | samtools view -b - \
  | samtools sort -@ 32 -o alignment/sample.sorted.bam
samtools index alignment/sample.sorted.bam
samtools flagstat alignment/sample.sorted.bam > alignment/sample.flagstat.txt
java -Xmx6g -jar $PICARD_JAR MarkDuplicates \
  I=alignment/sample.sorted.bam \
  O=alignment/sample.markdup.bam \
  M=alignment/sample.markdup.metrics.txt \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
```

### 3.4 BQSR（可选）
```bash
gatk BaseRecalibrator \
  -R ref-202509.fna \
  -I alignment/sample.markdup.bam \
  --known-sites known_sites.vcf.gz \
  -O gatk/sample.recal.table

gatk ApplyBQSR \
  -R ref-202509.fna \
  -I alignment/sample.markdup.bam \
  --bqsr-recal-file gatk/sample.recal.table \
  -O gatk/sample.bqsr.bam
```

### 3.5 变异检测与过滤
```bash
gatk HaplotypeCaller -R ref-202509.fna -I gatk/sample.bqsr.bam -O variants/sample.g.vcf.gz -ERC GVCF
gatk GenotypeGVCFs -R ref-202509.fna -V variants/sample.g.vcf.gz -O variants/sample.raw.vcf.gz
gatk VariantFiltration \
  -R ref-202509.fna \
  -V variants/sample.raw.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name BASIC_SNP_FILTER \
  -O variants/sample.tmp.filtered.vcf.gz
bcftools view -f PASS -Oz -o variants/sample.filtered.vcf.gz variants/sample.tmp.filtered.vcf.gz
vcftools --gzvcf variants/sample.filtered.vcf.gz --site-mean-depth --out variants/sample
vcftools --gzvcf variants/sample.filtered.vcf.gz --site-quality --out variants/sample
```

## 4. 可选扩展

- **多样本联合分析**：多个 `gVCF` 可通过 `gatk CombineGVCFs` 或 `GenomicsDBImport` 联合，再 `GenotypeGVCFs`。
- **结构变异检测**：可额外引入 `Manta`, `Delly` 等工具。
- **覆盖度评估**：`mosdepth`, `samtools depth` 可用于绘制覆盖度曲线。
- **注释**：建议使用 `ANNOVAR`, `SnpEff` 等对过滤后的 `VCF` 进行功能注释。

## 5. 常见问题

1. **参考基因组缺失索引**：脚本会自动尝试生成，但需确保对输出目录有写权限。
2. **BQSR 已知位点缺失**：可跳过 `-k` 参数，此时流程直接使用重复标记后的 BAM 进行调用。
3. **PICARD_JAR 路径**：若未指定 `-p`，脚本会尝试使用系统命令 `picard`。确保 `picard` 在 `PATH` 中或显式传入 jar 路径。
4. **T7 平台特有标签**：若需要在比对时添加 read group，可在 `bwa mem` 前添加 `-R '@RG\tID:...\tSM:...\tPL:BGISEQ'` 参数。

## 6. 结果解释

- `*.filtered.vcf.gz`：最终高质量 SNP/Indel 集合。
- `*.flagstat.txt`：比对统计，可检查 mapping rate、proper pair 比例等。
- `*.markdup.metrics.txt`：重复率统计，帮助评估建库质量。
- `vcftools` 输出：深度（`.ldepth.mean`）与质量（`.lqual`）统计，便于下游分析。

如需进一步自动化，可在现有脚本基础上改写为 `Snakemake` 或 `Nextflow` 流程。



---

# Auto-Fish-SNP-Calling: 自动化 SNP 变异检测流程

本项目详细介绍了一套用于高通量测序数据（兼容 Illumina 及华大 T7 平台）的自动化 SNP/Indel 变异检测流程。该流程通过两个核心脚本，实现了从环境一键部署到全自动分析的完整解决方案，覆盖了从原始 FASTQ 数据到高质量 VCF 文件的所有关键步骤。

-   `1_Set_Env.sh`: 自动化环境部署脚本，用于快速搭建包含所有必需软件的 Conda 环境。
-   `2_Main_Workflow.sh`: 核心分析流程脚本，支持单样本和批量处理模式，并提供丰富的参数化选项进行流程调优。

## 1. 环境准备与一键部署

在开始分析之前，您需要一个包含了所有必需生物信息学软件的运行环境。我们提供了 `1_Set_Env.sh` 脚本来自动完成这一过程。

**先决条件:**
*   操作系统：Linux
*   已安装 Conda 或 Miniconda，并已将其加入系统 `PATH`。
*   已安装 `wget` 和 `unzip` 工具。
*   稳定的网络连接。

**部署步骤:**

1.  **下载并执行部署脚本**
    为脚本添加执行权限并运行：
    ```bash
    chmod +x 1_Set_Env.sh
    ./1_Set_Env.sh
    ```
    该脚本将自动执行以下操作：
    *   检查系统依赖（Conda, wget, unzip）。
    *   配置 Conda 软件源（bioconda, conda-forge）。
    *   创建两个独立的 Conda 环境：
        *   `snp_call`: 主分析环境，包含 BWA, Samtools, Picard, GATK 依赖等。
        *   `fastqc`: 独立的质控环境，包含 FastQC。
    *   下载并安装 GATK（默认版本 4.6.2.0）。
    *   自动将 GATK 路径写入 `~/.bashrc` 配置文件。

2.  **激活环境变量**
    GATK 的环境变量配置需要重新加载 shell 配置才能生效。请执行以下**二选一**的操作：
    *   关闭当前终端，并重新打开一个新的终端。
    *   在当前终端中运行命令：`source ~/.bashrc`

3.  **激活分析环境**
    所有分析都应在 `snp_call` 环境中进行。在开始运行主流程前，请务必激活它：
    ```bash
    conda activate snp_call
    ```
    脚本内置了环境检查，若未在 `snp_call` 环境中运行，将会报错并提示。

## 2. 数据准备

1.  **原始测序数据**: 将成对的 `R1/R2` 原始测序文件（例如 `sampleA_R1.fastq.gz`, `sampleA_R2.fastq.gz`）准备好。批量处理时，建议将所有样本的 FASTQ 文件放入同一个输入目录。
2.  **参考基因组**: 准备 `FASTA` 格式的参考基因组文件（例如 `reference.fasta`）。脚本会自动为其创建 BWA, Samtools, GATK 所需的索引。
3.  **(可选) 已知变异位点**: 准备用于碱基质量再校正（BQSR）的已知变异位点 VCF 文件（例如 `known_sites.vcf.gz`）。这通常是来自 dbSNP 或 1000 Genomes Project 的资源。若不提供，流程将跳过 BQSR 步骤。

> **T7 数据兼容说明**：华大 T7 平台输出的 FASTQ 格式与 Illumina 完全兼容，可直接使用 `bwa`, `fastp` 等主流工具处理。

## 3. 快速开始：使用脚本一键运行

`2_Main_Workflow.sh` 脚本是流程的核心，它封装了所有分析步骤，并支持单样本和批量处理模式。

### 3.1 单样本处理模式

此模式适用于处理单个样本，需要明确指定输入输出文件。

```bash
# 示例：处理名为 sample01 的单个样本
bash 2_Main_Workflow.sh \
  -1 /path/to/data/sample01_R1.fastq.gz \
  -2 /path/to/data/sample01_R2.fastq.gz \
  -r /path/to/reference/ref.fna \
  -o /path/to/results/sample01 \
  -s sample01 \
  -k /path/to/known_sites/dbsnp.vcf.gz \
  -j 1 \
  -c 0.8
```

### 3.2 批量处理模式

此模式会自动发现指定目录下的所有成对 FASTQ 文件，并以多进程方式并行处理它们。

```bash
# 示例：并行处理 fastq_dir 目录下的所有样本，同时运行 4 个样本
bash 2_Main_Workflow.sh \
  -d /path/to/fastq_dir \
  -r /path/to/reference/ref.fna \
  -o /path/to/batch_results \
  -j 4 \
  -m 20 \
  -c 0.9
```

脚本将为 `fastq_dir` 中的每个样本（如 `sampleA`, `sampleB`）在输出目录 `/path/to/batch_results` 下创建独立的子目录（`sampleA/`, `sampleB/`）。

## 4. 完整参数参考

脚本提供了丰富的参数以满足不同场景下的分析需求。

| 分类 | 参数 | 描述 | 默认值 |
| :--- | :--- | :--- | :--- |
| **输入/输出** | `-r <file>` | **必填**，参考基因组 FASTA 文件。 | |
| | `-o <dir>` | **必填**，输出目录。单样本模式下为最终输出目录；批量模式下为总输出目录。 | |
| | `-1 <file>` | （单样本模式）Read1 FASTQ 文件（.gz 格式）。 | |
| | `-2 <file>` | （单样本模式）Read2 FASTQ 文件（.gz 格式）。 | |
| | `-s <string>` | （单样本模式）样本 ID，用于文件命名。 | |
| | `-d <dir>` | （批量模式）包含成对 FASTQ 文件的目录。 | |
| **资源管理** | `-j <int>` | 并行处理的样本组数（仅批量模式）。 | `1` |
| | `-m <int>` | 启动新任务所需的最小空闲内存（GB）（仅批量模式）。 | `10` |
| | `-c <float>` | CPU 使用率上限（0.1-1.0），用于计算流程可用的总线程数。 | `0.8` |
| **环境与工具**| `-q <name>` | 用于运行 FastQC 的独立 Conda 环境名称。 | （空，使用当前环境）|
| | `-p <file>` | Picard jar 文件路径。若为空，则尝试直接调用 `picard` 命令。 | `""` |
| | `-G <file>` | GATK 可执行文件路径。 | `gatk` |
| **流程调优** | `-k <file>` | 已知变异位点 VCF 文件，用于 BQSR（可多次使用）。 | （空） |
| | `-g <string>`| BWA Read Group 模板字符串，`{sample}` 会被替换为样本ID。 | `""` |
| | `-L <file>` | GATK：仅分析指定区域（BED 或 Interval List 文件）。 | （空） |
| | `--fastp-qual <int>` | fastp: 平均质量值阈值。 | `20` |
| | `--fastp-len <int>` | fastp: 过滤后保留的最小读长。 | `50` |
| | `--bwa-use-M` | BWA-MEM: 添加 `-M` 标记（推荐）。 | `true` (开启) |
| | `--ploidy <int>` | GATK: 物种倍性。 | `2` |
| | `-f` | 强制重跑所有步骤，忽略已完成的状态记录。 | `false` |
| **变异过滤** | `--filter-qd <float>` | 过滤 QD < value 的位点。 | `1.5` |
| | `--filter-fs <float>` | 过滤 FS > value 的位点。 | `60.0` |
| | `--filter-mq <float>` | 过滤 MQ < value 的位点。 | `30.0` |
| | `--filter-mqrs <float>` | 过滤 MQRankSum < value 的位点。 | `-12.5` |
| | `--filter-rprs <float>`| 过滤 ReadPosRankSum < value 的位点。 | `-8.0` |
| | `--filter-sor <float>` | 过滤 SOR > value 的位点。 | `3.0` |
| | `--filter-dp <int>` | 过滤 DP < value 的位点。 | `10` |

## 5. 流程步骤与输出结构

脚本会自动执行 GATK Best Practices 的核心步骤，所有中间结果和日志均在输出目录下分类保存。

**核心分析步骤:**
1.  **参考基因组索引**: 自动检查并创建 BWA, Samtools (.fai) 和 GATK (.dict) 索引。
2.  **原始数据质控**: 使用 `fastqc` 对原始 FASTQ 文件进行质量评估。
3.  **数据清洗**: 使用 `fastp` 进行接头去除、低质量过滤和长度过滤。
4.  **清洗后质控**: 再次使用 `fastqc` 评估清洗后数据的质量。
5.  **序列比对**: 使用 `bwa mem` 将 reads 比对到参考基因组，并使用 `samtools` 排序和索引，生成 BAM 文件。
6.  **标记重复**: 使用 `picard MarkDuplicates` 标记 PCR 重复序列。
7.  **碱基质量再校正 (BQSR)**: (可选) 如果提供了 `-k` 参数，则使用 `gatk BaseRecalibrator` 和 `ApplyBQSR` 修正测序系统误差。
8.  **变异检测**: 使用 `gatk HaplotypeCaller` 以 GVCF 模式进行变异检测。
9.  **基因分型**: 使用 `gatk GenotypeGVCFs` 将 GVCF 文件转换为 VCF 文件。
10. **变异过滤**: 使用 `gatk VariantFiltration` 对 VCF 进行硬过滤，并用 `bcftools` 提取通过过滤（PASS）的位点。
11. **统计分析**: 使用 `vcftools` 对最终的 VCF 文件进行深度和质量统计。

**输出目录结构 (`-o` 指定的目录):**
```
<output_dir>/
├── logs/                 # 主流程运行日志 (pipeline.log)
├── qc/
│   ├── raw/              # 原始 FASTQ 的 FastQC 报告
│   └── trimmed/          # 清洗后 FASTQ 的 FastQC 报告
├── trimmed/              # fastp 清洗后的 FASTQ 文件及 JSON/HTML 报告
├── alignment/            # 比对结果
│   ├── *.sorted.bam      # 排序后的 BAM 文件
│   ├── *.markdup.bam     # 标记重复后的 BAM 文件
│   ├── *.flagstat.txt    # 比对统计
│   └── *.metrics.txt     # 重复标记统计
├── gatk/                 # GATK 中间文件 (如 BQSR recal table)
└── variants/             # 变异检测结果
    ├── *.g.vcf.gz        # HaplotypeCaller 输出的 GVCF 文件
    ├── *.raw.vcf.gz      # GenotypeGVCFs 输出的原始 VCF
    ├── *.filtered.vcf.gz # 最终过滤后的高质量 VCF 文件
    └── *.ldepth.mean     # VCF 位点平均深度统计
    └── *.lqual           # VCF 位点质量统计
```

## 6. 常见问题 (FAQ)

1.  **启动脚本时提示 `conda: command not found` 或 `未检测到conda环境`**
    *   请确保 Conda/Miniconda 已正确安装并将其 `bin` 目录添加到了系统 `PATH`。
    *   在运行 `2_Main_Workflow.sh` 之前，务必使用 `conda activate snp_call` 命令激活正确的分析环境。

2.  **GATK 或 Picard 命令找不到**
    *   对于 GATK，请确认 `1_Set_Env.sh` 已成功运行，并且您已经通过 `source ~/.bashrc` 或重开终端的方式使环境变量生效。
    *   对于 Picard，如果 `picard` 命令不可用，请使用 `-p /path/to/picard.jar` 参数明确指定其 jar 文件路径。

3.  **BQSR 步骤是否必须？**
    *   不是。如果您没有高质量的已知变异位点数据库，可以不提供 `-k` 参数。流程会自动跳过 BQSR，直接使用标记重复后的 BAM 文件进行变异检测，这在处理非模式生物时是常见做法。

4.  **如何调整过滤参数？**
    *   GATK 的硬过滤阈值对结果影响显著。脚本提供了一系列 `--filter-*` 参数，您可以根据物种特性、测序深度和数据质量分布进行调整，以达到最佳的特异性和灵敏度平衡。建议先用默认值运行，然后根据 `VariantFiltration` 前后 VCF 的变化来优化参数。

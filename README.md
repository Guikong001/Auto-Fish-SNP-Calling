# Auto-Fish-SNP-Calling

本仓库提供一个适用于华大 T7/Illumina 测序平台的标准 SNP/Indel 变异检测流程示例，包括：

- `scripts/run_pipeline.sh`：一键式 Bash 脚本，可完成质控、清洗、比对、重复标记、可选 BQSR、变异检测与过滤。
- `docs/pipeline_overview.md`：详细的流程说明、手动执行命令以及常见问题解答。
- `requirements.txt`：所需生信软件列表（需自行安装）。

## 快速开始

1. 按 `requirements.txt` 安装所需软件，并确保命令可在当前环境调用。
2. 准备原始测序数据 `sample_R1.fastq.gz`、`sample_R2.fastq.gz` 和参考基因组 `ref-202509.fna`。
3. 执行示例命令：

```bash
bash scripts/run_pipeline.sh \
  -1 /path/to/sample_R1.fastq.gz \
  -2 /path/to/sample_R2.fastq.gz \
  -r /path/to/ref-202509.fna \
  -o results/sample01 \
  -s sample01 \
  -p /home/wslll/13_anti_disease_natural/0_biosoftware/picard.jar \
  -g '@RG\tID:lane1\tSM:sample01\tPL:BGISEQ'
```

如需使用 BQSR，请通过 `-k` 传入一个或多个已知变异位点 VCF；若需自定义 read group，可使用 `-g` 选项。

更多细节、手动执行流程、扩展建议请参阅 [docs/pipeline_overview.md](docs/pipeline_overview.md)。

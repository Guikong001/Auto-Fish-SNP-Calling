#!/usr/bin/env bash
#
# 自动化的T7/Illumina测序数据变异检测流程
#
set -euo pipefail

usage() {
  cat <<USAGE
用法: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -r <ref.fna> -o <out_dir> -s <sample_id> [选项]

必填参数:
  -1    Read1 FASTQ 文件 (gzip 压缩)
  -2    Read2 FASTQ 文件 (gzip 压缩)
  -r    参考基因组 fasta 文件 (例如 ref-202509.fna)
  -o    输出目录
  -s    样本 ID (用于文件命名)

可选参数:
  -k    已知变异位点 (用于BQSR，可多次指定)
  -t    线程数 (默认: 16)
  -p    PICARD_JAR 路径 (默认读取环境变量 PICARD_JAR 或系统中的 picard 命令)
  -g    BWA read group 信息 (示例: '@RG\tID:lane1\tSM:sample\tPL:BGISEQ')
  -h    显示此帮助信息

示例:
  $0 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -r ref-202509.fna \
     -o results/sampleA -s sampleA -k known_sites.vcf.gz -t 32
USAGE
}

RAW_R1=""
RAW_R2=""
REF=""
OUTDIR=""
SAMPLE=""
THREADS=16
KNOWN_SITES=()
PICARD_PATH="${PICARD_JAR:-}"
READ_GROUP=""

while getopts ":1:2:r:o:s:k:t:p:g:h" opt; do
  case $opt in
    1) RAW_R1=$OPTARG ;;
    2) RAW_R2=$OPTARG ;;
    r) REF=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;
    k) KNOWN_SITES+=("$OPTARG") ;;
    t) THREADS=$OPTARG ;;
    p) PICARD_PATH=$OPTARG ;;
    g) READ_GROUP=$OPTARG ;;
    h) usage; exit 0 ;;
    :) echo "参数 -$OPTARG 需要一个值" >&2; usage; exit 1 ;;
    \?) echo "未知参数: -$OPTARG" >&2; usage; exit 1 ;;
  esac
done

if [[ -z $RAW_R1 || -z $RAW_R2 || -z $REF || -z $OUTDIR || -z $SAMPLE ]]; then
  echo "缺少必填参数" >&2
  usage
  exit 1
fi

for f in "$RAW_R1" "$RAW_R2" "$REF"; do
  if [[ ! -f $f ]]; then
    echo "文件不存在: $f" >&2
    exit 1
  fi
done

# 检查必要工具
REQUIRED_CMDS=(fastqc fastp bwa samtools gatk bcftools vcftools)
for cmd in "${REQUIRED_CMDS[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "缺少依赖命令: $cmd" >&2
    exit 1
  fi
done

if [[ -n $PICARD_PATH && ! -f $PICARD_PATH ]]; then
  echo "PICARD_JAR 文件不存在: $PICARD_PATH" >&2
  exit 1
fi

mkdir -p "$OUTDIR"
LOG_DIR="$OUTDIR/logs"
RAW_QC_DIR="$OUTDIR/qc/raw"
TRIM_QC_DIR="$OUTDIR/qc/trimmed"
TRIM_DIR="$OUTDIR/trimmed"
ALIGN_DIR="$OUTDIR/alignment"
GATK_DIR="$OUTDIR/gatk"
VARIANT_DIR="$OUTDIR/variants"
REPORT_DIR="$OUTDIR/report"

mkdir -p "$LOG_DIR" "$RAW_QC_DIR" "$TRIM_QC_DIR" "$TRIM_DIR" \
         "$ALIGN_DIR" "$GATK_DIR" "$VARIANT_DIR" "$REPORT_DIR"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

run_picard() {
  local tool=$1; shift
  if [[ -n $PICARD_PATH ]]; then
    java -Xmx6g -jar "$PICARD_PATH" "$tool" "$@"
  else
    picard "$tool" "$@"
  fi
}

# 1. 原始数据质控
log "运行 FastQC (raw reads)"
fastqc -t "$THREADS" -o "$RAW_QC_DIR" "$RAW_R1" "$RAW_R2"

# 2. fastp 清洗
TRIMMED_R1="$TRIM_DIR/${SAMPLE}_R1.trimmed.fastq.gz"
TRIMMED_R2="$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz"
FASTP_JSON="$TRIM_DIR/${SAMPLE}.fastp.json"
FASTP_HTML="$TRIM_DIR/${SAMPLE}.fastp.html"
log "运行 fastp 进行质控/剪切"
fastp \
  --thread "$THREADS" \
  --detect_adapter_for_pe \
  --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 50 \
  -i "$RAW_R1" -I "$RAW_R2" \
  -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
  --json "$FASTP_JSON" --html "$FASTP_HTML" \
  --report_title "${SAMPLE} fastp report"

log "运行 FastQC (trimmed reads)"
fastqc -t "$THREADS" -o "$TRIM_QC_DIR" "$TRIMMED_R1" "$TRIMMED_R2"

# 3. 整理质控报告
if command -v multiqc >/dev/null 2>&1; then
  log "生成 MultiQC 汇总报告"
  multiqc "$OUTDIR" -o "$REPORT_DIR"
fi

# 4. 准备参考基因组索引
log "检查并生成参考基因组索引"
if [[ ! -f "${REF}.bwt" && ! -f "${REF}.0123" ]]; then
  bwa index "$REF"
fi
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi
REF_DICT="${REF%.*}.dict"
if [[ ! -f $REF_DICT ]]; then
  gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"
fi

# 5. 比对并排序
SORTED_BAM="$ALIGN_DIR/${SAMPLE}.sorted.bam"
log "运行 BWA-MEM 比对并排序"
if [[ -n $READ_GROUP ]]; then
  bwa mem -t "$THREADS" -R "$READ_GROUP" "$REF" "$TRIMMED_R1" "$TRIMMED_R2" \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM"
else
  bwa mem -t "$THREADS" "$REF" "$TRIMMED_R1" "$TRIMMED_R2" \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM"
fi

samtools index "$SORTED_BAM"

samtools flagstat "$SORTED_BAM" > "$ALIGN_DIR/${SAMPLE}.flagstat.txt"

# 6. 标记重复
MARKED_BAM="$ALIGN_DIR/${SAMPLE}.markdup.bam"
METRICS_FILE="$ALIGN_DIR/${SAMPLE}.markdup.metrics.txt"
log "运行 Picard MarkDuplicates"
run_picard MarkDuplicates \
  I="$SORTED_BAM" \
  O="$MARKED_BAM" \
  M="$METRICS_FILE" \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT

BAM_FOR_CALLING="$MARKED_BAM"

# 7. (可选) BQSR
if [[ ${#KNOWN_SITES[@]} -gt 0 ]]; then
  RECAL_TABLE="$GATK_DIR/${SAMPLE}.recal.table"
  BQSR_BAM="$GATK_DIR/${SAMPLE}.bqsr.bam"
  log "运行 GATK BaseRecalibrator"
  gatk BaseRecalibrator \
    -R "$REF" \
    -I "$MARKED_BAM" \
    $(printf ' --known-sites %q' "${KNOWN_SITES[@]}") \
    -O "$RECAL_TABLE"

  log "应用 BQSR 校正"
  gatk ApplyBQSR \
    -R "$REF" \
    -I "$MARKED_BAM" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$BQSR_BAM"

  samtools index "$BQSR_BAM"
  BAM_FOR_CALLING="$BQSR_BAM"
fi

# 8. 变异检测
GVCF="$VARIANT_DIR/${SAMPLE}.g.vcf.gz"
RAW_VCF="$VARIANT_DIR/${SAMPLE}.raw.vcf.gz"
FILTERED_VCF="$VARIANT_DIR/${SAMPLE}.filtered.vcf.gz"

log "运行 GATK HaplotypeCaller (GVCF 模式)"
gatk HaplotypeCaller \
  -R "$REF" \
  -I "$BAM_FOR_CALLING" \
  -O "$GVCF" \
  -ERC GVCF

log "将 GVCF 转换为常规 VCF"
gatk GenotypeGVCFs \
  -R "$REF" \
  -V "$GVCF" \
  -O "$RAW_VCF"

log "执行变异过滤"
gatk VariantFiltration \
  -R "$REF" \
  -V "$RAW_VCF" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "BASIC_SNP_FILTER" \
  -O "$VARIANT_DIR/${SAMPLE}.tmp.filtered.vcf.gz"

bcftools view -f PASS -Oz -o "$FILTERED_VCF" "$VARIANT_DIR/${SAMPLE}.tmp.filtered.vcf.gz"
if command -v tabix >/dev/null 2>&1; then
  tabix -p vcf "$FILTERED_VCF"
fi
rm -f "$VARIANT_DIR/${SAMPLE}.tmp.filtered.vcf.gz" "$VARIANT_DIR/${SAMPLE}.tmp.filtered.vcf.gz.tbi" 2>/dev/null || true

log "生成 VCF 统计信息"
vcftools --gzvcf "$FILTERED_VCF" --site-mean-depth --out "$VARIANT_DIR/${SAMPLE}"
vcftools --gzvcf "$FILTERED_VCF" --site-quality --out "$VARIANT_DIR/${SAMPLE}"

log "流程完成。最终结果:"
log "  变异结果: $FILTERED_VCF"
log "  统计信息: $VARIANT_DIR/${SAMPLE}.site-mean-depth"
log "  统计信息: $VARIANT_DIR/${SAMPLE}.site-quality"

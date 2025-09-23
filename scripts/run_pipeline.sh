#!/usr/bin/env bash
#
# 自动化的T7/Illumina测序数据变异检测流程
#
set -euo pipefail

usage() {
  cat <<USAGE
用法:
  单样本: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -r <ref.fna> -o <out_dir> -s <sample_id> [选项]
  批量处理目录: $0 -d <fastq_dir> -r <ref.fna> -o <out_dir> [选项]

必填参数:
  -r    参考基因组 fasta 文件 (例如 ref-202509.fna)
  -o    输出目录 (批量处理时会在此目录下生成按样本划分的子目录)

输入数据（二选一）:
  -1    Read1 FASTQ 文件 (gzip 压缩)
  -2    Read2 FASTQ 文件 (gzip 压缩)
  -s    样本 ID (用于单样本模式的文件命名)
  -d    成对 FASTQ 所在目录，文件命名形如 586_R1.fq.gz / 586_R2.fq.gz

可选参数:
  -k    已知变异位点 (用于BQSR，可多次指定)
  -t    线程数 (默认: 16)
  -p    PICARD_JAR 路径 (默认读取环境变量 PICARD_JAR 或系统中的 picard 命令)
  -g    BWA read group 模板 (支持占位符 {sample})
  -G    GATK 可执行文件路径 (默认使用系统中的 gatk 命令)
  -h    显示此帮助信息

示例:
  单样本: $0 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -r ref-202509.fna \
            -o results/sampleA -s sampleA -k known_sites.vcf.gz -t 32 -G /home/wslll/GATK-4.6.2/gatk-4.6.2.0/gatk
  批量:   $0 -d fastq_dir -r ref-202509.fna -o results -t 32
USAGE
}

RAW_R1=""
RAW_R2=""
REF=""
OUTDIR=""
SAMPLE=""
DATA_DIR=""
THREADS=16
KNOWN_SITES=()
PICARD_PATH="${PICARD_JAR:-}"
READ_GROUP_TEMPLATE=""
GATK_CMD="${GATK_PATH:-gatk}"

while getopts ":1:2:r:o:s:d:k:t:p:g:G:h" opt; do
  case $opt in
    1) RAW_R1=$OPTARG ;;
    2) RAW_R2=$OPTARG ;;
    r) REF=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;
    d) DATA_DIR=$OPTARG ;;
    k) KNOWN_SITES+=("$OPTARG") ;;
    t) THREADS=$OPTARG ;;
    p) PICARD_PATH=$OPTARG ;;
    g) READ_GROUP_TEMPLATE=$OPTARG ;;
    G) GATK_CMD=$OPTARG ;;
    h) usage; exit 0 ;;
    :) echo "参数 -$OPTARG 需要一个值" >&2; usage; exit 1 ;;
    \?) echo "未知参数: -$OPTARG" >&2; usage; exit 1 ;;
  esac
done

if [[ -z $REF || -z $OUTDIR ]]; then
  echo "缺少必填参数" >&2
  usage
  exit 1
fi

if [[ -n $DATA_DIR ]]; then
  if [[ -n $RAW_R1 || -n $RAW_R2 || -n $SAMPLE ]]; then
    echo "批量模式下不需要 -1/-2/-s" >&2
    exit 1
  fi
  if [[ ! -d $DATA_DIR ]]; then
    echo "目录不存在: $DATA_DIR" >&2
    exit 1
  fi
else
  if [[ -z $RAW_R1 || -z $RAW_R2 || -z $SAMPLE ]]; then
    echo "缺少必填参数" >&2
    usage
    exit 1
  fi
  for f in "$RAW_R1" "$RAW_R2"; do
    if [[ ! -f $f ]]; then
      echo "文件不存在: $f" >&2
      exit 1
    fi
  done
fi

if [[ ! -f $REF ]]; then
  echo "参考基因组不存在: $REF" >&2
  exit 1
fi

# 检查必要工具
REQUIRED_CMDS=(fastqc fastp bwa samtools bcftools vcftools)
for cmd in "${REQUIRED_CMDS[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "缺少依赖命令: $cmd" >&2
    exit 1
  fi
done

if [[ $GATK_CMD == gatk ]]; then
  if ! command -v gatk >/dev/null 2>&1; then
    echo "缺少依赖命令: gatk" >&2
    exit 1
  fi
else
  if [[ ! -x $GATK_CMD ]]; then
    echo "GATK 可执行文件不可用: $GATK_CMD" >&2
    exit 1
  fi
fi

if [[ -n $PICARD_PATH && ! -f $PICARD_PATH ]]; then
  echo "PICARD_JAR 文件不存在: $PICARD_PATH" >&2
  exit 1
fi

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

run_gatk() {
  "$GATK_CMD" "$@"
}

setup_dirs() {
  local base_dir=$1
  mkdir -p "$base_dir"
  LOG_DIR="$base_dir/logs"
  RAW_QC_DIR="$base_dir/qc/raw"
  TRIM_QC_DIR="$base_dir/qc/trimmed"
  TRIM_DIR="$base_dir/trimmed"
  ALIGN_DIR="$base_dir/alignment"
  GATK_DIR="$base_dir/gatk"
  VARIANT_DIR="$base_dir/variants"
  REPORT_DIR="$base_dir/report"
  mkdir -p "$LOG_DIR" "$RAW_QC_DIR" "$TRIM_QC_DIR" "$TRIM_DIR" \
           "$ALIGN_DIR" "$GATK_DIR" "$VARIANT_DIR" "$REPORT_DIR"
}

resolve_read_group() {
  local sample_id=$1
  local rg="$READ_GROUP_TEMPLATE"
  if [[ -z $rg ]]; then
    rg="@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA"
  else
    rg=${rg//\{sample\}/$sample_id}
  fi
  printf '%s' "$rg"
}

process_sample() {
  local sample=$1
  local raw_r1=$2
  local raw_r2=$3
  local sample_out=$4

  log "==== 开始处理样本 ${sample} ===="
  setup_dirs "$sample_out"

  for f in "$raw_r1" "$raw_r2"; do
    if [[ ! -f $f ]]; then
      echo "文件不存在: $f" >&2
      return 1
    fi
  done

# 1. 原始数据质控
  log "运行 FastQC (raw reads)"
  fastqc -t "$THREADS" -o "$RAW_QC_DIR" "$raw_r1" "$raw_r2"

# 2. fastp 清洗
  TRIMMED_R1="$TRIM_DIR/${sample}_R1.trimmed.fastq.gz"
  TRIMMED_R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
  FASTP_JSON="$TRIM_DIR/${sample}.fastp.json"
  FASTP_HTML="$TRIM_DIR/${sample}.fastp.html"
  log "运行 fastp 进行质控/剪切"
  fastp \
    --thread "$THREADS" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 \
    -i "$raw_r1" -I "$raw_r2" \
    -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" \
    --report_title "${sample} fastp report"

  log "运行 FastQC (trimmed reads)"
  fastqc -t "$THREADS" -o "$TRIM_QC_DIR" "$TRIMMED_R1" "$TRIMMED_R2"

# 3. 整理质控报告
  if command -v multiqc >/dev/null 2>&1; then
    log "生成 MultiQC 汇总报告"
    multiqc "$sample_out" -o "$REPORT_DIR"
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
    run_gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"
  fi

# 5. 比对并排序
  SORTED_BAM="$ALIGN_DIR/${sample}.sorted.bam"
  log "运行 BWA-MEM 比对并排序"
  local read_group
  read_group=$(resolve_read_group "$sample")
  bwa mem -t "$THREADS" -R "$read_group" "$REF" "$TRIMMED_R1" "$TRIMMED_R2" \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM"

  samtools index "$SORTED_BAM"

  samtools flagstat "$SORTED_BAM" > "$ALIGN_DIR/${sample}.flagstat.txt"

# 6. 标记重复
  MARKED_BAM="$ALIGN_DIR/${sample}.markdup.bam"
  METRICS_FILE="$ALIGN_DIR/${sample}.markdup.metrics.txt"
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
    RECAL_TABLE="$GATK_DIR/${sample}.recal.table"
    BQSR_BAM="$GATK_DIR/${sample}.bqsr.bam"
    log "运行 GATK BaseRecalibrator"
    run_gatk BaseRecalibrator \
      -R "$REF" \
      -I "$MARKED_BAM" \
      $(printf ' --known-sites %q' "${KNOWN_SITES[@]}") \
      -O "$RECAL_TABLE"

    log "应用 BQSR 校正"
    run_gatk ApplyBQSR \
      -R "$REF" \
      -I "$MARKED_BAM" \
      --bqsr-recal-file "$RECAL_TABLE" \
      -O "$BQSR_BAM"

    samtools index "$BQSR_BAM"
    BAM_FOR_CALLING="$BQSR_BAM"
  fi

# 8. 变异检测
  GVCF="$VARIANT_DIR/${sample}.g.vcf.gz"
  RAW_VCF="$VARIANT_DIR/${sample}.raw.vcf.gz"
  FILTERED_VCF="$VARIANT_DIR/${sample}.filtered.vcf.gz"

  log "运行 GATK HaplotypeCaller (GVCF 模式)"
  run_gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM_FOR_CALLING" \
    -O "$GVCF" \
    -ERC GVCF

  log "将 GVCF 转换为常规 VCF"
  run_gatk GenotypeGVCFs \
    -R "$REF" \
    -V "$GVCF" \
    -O "$RAW_VCF"

  log "执行变异过滤"
  run_gatk VariantFiltration \
    -R "$REF" \
    -V "$RAW_VCF" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "BASIC_SNP_FILTER" \
    -O "$VARIANT_DIR/${sample}.tmp.filtered.vcf.gz"

  bcftools view -f PASS -Oz -o "$FILTERED_VCF" "$VARIANT_DIR/${sample}.tmp.filtered.vcf.gz"
  if command -v tabix >/dev/null 2>&1; then
    tabix -p vcf "$FILTERED_VCF"
  fi
  rm -f "$VARIANT_DIR/${sample}.tmp.filtered.vcf.gz" "$VARIANT_DIR/${sample}.tmp.filtered.vcf.gz.tbi" 2>/dev/null || true

  log "生成 VCF 统计信息"
  vcftools --gzvcf "$FILTERED_VCF" --site-mean-depth --out "$VARIANT_DIR/${sample}"
  vcftools --gzvcf "$FILTERED_VCF" --site-quality --out "$VARIANT_DIR/${sample}"

  log "流程完成。最终结果:"
  log "  变异结果: $FILTERED_VCF"
  log "  统计信息: $VARIANT_DIR/${sample}.site-mean-depth"
  log "  统计信息: $VARIANT_DIR/${sample}.site-quality"
  log "==== 样本 ${sample} 处理完成 ===="
}

detect_pairs() {
  local dir=$1
  shopt -s nullglob
  local files=("$dir"/*_R1.fastq.gz "$dir"/*_R1.fq.gz "$dir"/*_1.fastq.gz "$dir"/*_1.fq.gz)
  shopt -u nullglob
  printf '%s\n' "${files[@]}"
}

find_r2_for() {
  local r1=$1
  local r2=""
  if [[ -f ${r1/_R1.fastq.gz/_R2.fastq.gz} ]]; then
    r2=${r1/_R1.fastq.gz/_R2.fastq.gz}
  elif [[ -f ${r1/_R1.fq.gz/_R2.fq.gz} ]]; then
    r2=${r1/_R1.fq.gz/_R2.fq.gz}
  elif [[ -f ${r1/_1.fastq.gz/_2.fastq.gz} ]]; then
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
  elif [[ -f ${r1/_1.fq.gz/_2.fq.gz} ]]; then
    r2=${r1/_1.fq.gz/_2.fq.gz}
  fi
  printf '%s' "$r2"
}

derive_sample_id() {
  local r1_basename=$1
  local sample=${r1_basename%%_R1.fastq.gz}
  if [[ $sample == "$r1_basename" ]]; then
    sample=${r1_basename%%_R1.fq.gz}
  fi
  if [[ $sample == "$r1_basename" ]]; then
    sample=${r1_basename%%_1.fastq.gz}
  fi
  if [[ $sample == "$r1_basename" ]]; then
    sample=${r1_basename%%_1.fq.gz}
  fi
  if [[ $sample == "$r1_basename" ]]; then
    sample=${sample%%.fastq.gz}
  fi
  if [[ $sample == "$r1_basename" ]]; then
    sample=${sample%%.fq.gz}
  fi
  printf '%s' "$sample"
}

if [[ -n $DATA_DIR ]]; then
  mkdir -p "$OUTDIR"
  mapfile -t R1_FILES < <(detect_pairs "$DATA_DIR")
  if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "未在目录 $DATA_DIR 中检测到 *_R1.fq.gz / *_R1.fastq.gz 文件" >&2
    exit 1
  fi
  for r1 in "${R1_FILES[@]}"; do
    r2=$(find_r2_for "$r1")
    if [[ -z $r2 ]]; then
      log "跳过，未找到匹配的 R2: $r1"
      continue
    fi
    sample=$(derive_sample_id "$(basename "$r1")")
    if [[ -z $sample ]]; then
      log "跳过，无法解析样本 ID: $r1"
      continue
    fi
    process_sample "$sample" "$r1" "$r2" "$OUTDIR/$sample"
  done
else
  process_sample "$SAMPLE" "$RAW_R1" "$RAW_R2" "$OUTDIR"
fi

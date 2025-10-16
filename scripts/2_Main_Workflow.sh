#!/usr/bin/env bash
#
# 自动化的T7/Illumina测序数据变异检测流程
# V4.1 - 全面参数化：fastp, BWA, GATK HaplotypeCaller 参数均可配置
#
set -euo pipefail

# --- 函数定义区域 ---

check_conda_environment() {
  # 检查是否在conda环境中
  if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
    echo "错误: 未检测到conda环境。请先激活conda环境。" >&2
    echo "建议运行: conda activate snp_call" >&2
    exit 1
  fi
  
  # 检查当前环境是否为snp_call
  if [[ "$CONDA_DEFAULT_ENV" != "snp_call" ]]; then
    echo "错误: 当前conda环境为 '$CONDA_DEFAULT_ENV'，但此脚本需要在 'snp_call' 环境中运行。" >&2
    echo "请运行以下命令切换到正确的环境:" >&2
    echo "  conda activate snp_call" >&2
    exit 1
  fi
  
  log "✓ conda环境检查通过: 当前环境为 '$CONDA_DEFAULT_ENV'"
}

usage() {
  cat <<USAGE
用法:
  单样本: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -r <ref.fna> -o <out_dir> -s <sample_id> [选项]
  批量处理: $0 -d <fastq_dir> -r <ref.fna> -o <out_dir> [选项]

必填参数:
  -r    参考基因组 fasta 文件
  -o    输出目录

输入数据（二选一）:
  -1    Read1 FASTQ 文件 (gzip 压缩)
  -2    Read2 FASTQ 文件 (gzip 压缩)
  -s    样本 ID (用于单样本模式的文件命名)
  -d    成对 FASTQ 所在目录

资源管理选项:
  -j    并行处理的样本组数 (默认: 1)
  -m    启动新任务所需最小空闲内存(GB) (默认: 10)
  -c    CPU使用率 (0.1-1.0), 用于计算流程可用的总线程数 (默认: 0.8)

环境管理选项:
  -q    用于运行 FastQC 的独立 Conda 环境名称 (例如: fastq)

##### 优化点 1: 在帮助信息中添加所有新的可调参数 #####
流程细节调优选项:
  --fastp-qual <int>      fastp: 平均质量值阈值 (默认: 20)
  --fastp-len <int>       fastp: 最小读长要求 (默认: 50)
  --bwa-use-M             BWA-MEM: 添加 -M 标记 (推荐, 默认开启)
  --ploidy <int>          GATK: 物种倍性 (默认: 2)
  -L <file>               GATK: 仅分析指定区域 (BED或Interval List文件)

过滤选项 (GATK VariantFiltration):
  --filter-qd <float>   过滤 QD < value 的位点 (默认: 1.5)
  --filter-fs <float>   过滤 FS > value 的位点 (默认: 60.0)
  --filter-mq <float>   过滤 MQ < value 的位点 (默认: 30.0)
  --filter-mqrs <float> 过滤 MQRankSum < value 的位点 (默认: -12.5)
  --filter-rprs <float> 过滤 ReadPosRankSum < value 的位点 (默认: -8.0)
  --filter-sor <float>  过滤 SOR > value 的位点 (默认: 3.0)
  --filter-dp <int>     过滤 DP < value 的位点 (默认: 10)

其他可选参数:
  -k    已知变异位点 (用于BQSR)
  -p    PICARD_JAR 路径
  -g    BWA read group 模板
  -G    GATK 可执行文件路径
  -f    强制重跑
  -h    显示此帮助信息
USAGE
}

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_memory_and_wait() {
    local required_gb=$1
    if (( required_gb <= 0 )); then return 0; fi
    local required_mb=$((required_gb * 1024))
    
    while true; do
        local available_mem_kb
        available_mem_kb=$(grep 'MemAvailable:' /proc/meminfo | awk '{print $2}')

        if [[ -z "$available_mem_kb" ]]; then
            log "警告: 'MemAvailable' 在 /proc/meminfo 中未找到，将使用 'Free + Buffers + Cached' 作为估算值。"
            local free_kb=$(grep 'MemFree:' /proc/meminfo | awk '{print $2}')
            local buffers_kb=$(grep 'Buffers:' /proc/meminfo | awk '{print $2}')
            local cached_kb=$(grep 'Cached:' /proc/meminfo | awk '{print $2}')
            available_mem_kb=$((free_kb + buffers_kb + cached_kb))
        fi

        local available_mem_mb=$((available_mem_kb / 1024))

        if (( available_mem_mb >= required_mb )); then
            log "可用内存 (${available_mem_mb}MB) 充足, 继续..."
            break
        else
            log "可用内存 (${available_mem_mb}MB) 低于阈值 (${required_mb}MB), 等待60秒后重试..."
            sleep 60
        fi
    done
}


run_or_skip() {
  local stage_name=$1; local target_file=$2; local command_str=$3; local status_file=$4; local log_prefix=$5
  local success_marker="STAGE_COMPLETE: ${stage_name}"
  if [[ $FORCE_RERUN == false ]] && grep -qFx "$success_marker" "$status_file" 2>/dev/null && [[ -s "$target_file" ]]; then
    log "${log_prefix} 跳过步骤 [${stage_name}] - 已完成"
    return 0
  fi
  log "${log_prefix} ==> 开始步骤 [${stage_name}]"
  if eval "$command_str"; then
    log "${log_prefix} ==> 步骤 [${stage_name}] 成功完成"
    echo "$success_marker" >> "$status_file"
  else
    log "${log_prefix} !!! 错误: 步骤 [${stage_name}] 失败"
    return 1
  fi
}

run_picard() {
  local tool=$1; shift
  if [[ -n $PICARD_PATH ]]; then
    java -Xmx8g -jar "$PICARD_PATH" "$tool" "$@"
  else
    picard "$tool" "$@"
  fi
}

run_gatk() {
  "$GATK_CMD" "$@"
}

setup_dirs() {
  mkdir -p "$1"/{logs,qc/{raw,trimmed},trimmed,alignment,gatk,variants,report}
}

resolve_read_group() {
  local sample_id=$1
  local rg="$READ_GROUP_TEMPLATE"
  if [[ -z $rg ]]; then
    rg="@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
  else
    rg=${rg//\{sample\}/$sample_id}
  fi
  printf '%s' "$rg"
}

process_sample() {
  local sample=$1; local raw_r1=$2; local raw_r2=$3; local sample_out=$4; local threads_to_use=$5
  local LOG_PREFIX="[样本: ${sample}]"
  log "${LOG_PREFIX} PID:$$ ==== 开始处理 ===="
  setup_dirs "$sample_out"
  local STATUS_FILE="$sample_out/.pipeline_status.log"
  if [[ $FORCE_RERUN == true && -f "$STATUS_FILE" ]]; then
    log "${LOG_PREFIX} 强制重跑，删除旧的状态记录"
    rm "$STATUS_FILE"
  fi
  touch "$STATUS_FILE"

  local qc_cmd_prefix=""
  if [[ -n "$QC_ENV_NAME" ]]; then
    qc_cmd_prefix="conda run -n \"$QC_ENV_NAME\""
  fi

  # 1. 原始数据质控
  local raw_base; raw_base=$(basename "${raw_r1}"); raw_base=${raw_base%.f*q.gz}
  local QC_RAW_DONE_MARKER="$sample_out/qc/raw/${raw_base}_fastqc.html"
  local cmd_qc_raw="${qc_cmd_prefix} fastqc -t \"$threads_to_use\" -o \"$sample_out/qc/raw\" \"$raw_r1\" \"$raw_r2\""
  run_or_skip "QC_RAW" "$QC_RAW_DONE_MARKER" "$cmd_qc_raw" "$STATUS_FILE" "$LOG_PREFIX"

  ##### 优化点 5: 在命令中动态使用 fastp 参数 #####
  # 2. fastp 清洗
  local TRIM_DIR="$sample_out/trimmed"
  local TRIMMED_R1="$TRIM_DIR/${sample}_R1.trimmed.fastq.gz"
  local TRIMMED_R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
  local cmd_fastp="fastp --thread \"$threads_to_use\" --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality \"${FASTP_QUAL}\" --length_required \"${FASTP_LEN}\" -i \"$raw_r1\" -I \"$raw_r2\" -o \"$TRIMMED_R1\" -O \"$TRIMMED_R2\" --json \"$TRIM_DIR/${sample}.fastp.json\" --html \"$TRIM_DIR/${sample}.fastp.html\" --report_title \"${sample} fastp report\""
  run_or_skip "TRIMMING" "$TRIMMED_R1" "$cmd_fastp" "$STATUS_FILE" "$LOG_PREFIX"

  # 3. 清洗后数据质控
  local trim_base; trim_base=$(basename "${TRIMMED_R1}"); trim_base=${trim_base%.f*q.gz}
  local QC_TRIM_DONE_MARKER="$sample_out/qc/trimmed/${trim_base}_fastqc.html"
  local cmd_qc_trim="${qc_cmd_prefix} fastqc -t \"$threads_to_use\" -o \"$sample_out/qc/trimmed\" \"$TRIMMED_R1\" \"$TRIMMED_R2\""
  run_or_skip "QC_TRIMMED" "$QC_TRIM_DONE_MARKER" "$cmd_qc_trim" "$STATUS_FILE" "$LOG_PREFIX"

  ##### 优化点 6: 在命令中动态使用 BWA-MEM 参数 #####
  # 5. 比对并排序
  local ALIGN_DIR="$sample_out/alignment"
  local SORTED_BAM="$ALIGN_DIR/${sample}.sorted.bam"
  local read_group; read_group=$(resolve_read_group "$sample")
  local bwa_opts=""
  if [[ "$BWA_USE_M" == true ]]; then bwa_opts="-M"; fi
  local cmd_align="bwa mem ${bwa_opts} -t \"$threads_to_use\" -R \"$read_group\" \"$REF\" \"$TRIMMED_R1\" \"$TRIMMED_R2\" | samtools view -b - | samtools sort -@ \"$threads_to_use\" -o \"$SORTED_BAM\" && samtools index \"$SORTED_BAM\" && samtools flagstat \"$SORTED_BAM\" > \"$ALIGN_DIR/${sample}.flagstat.txt\""
  run_or_skip "ALIGNMENT" "$SORTED_BAM" "$cmd_align" "$STATUS_FILE" "$LOG_PREFIX"
  
  # ... 后续步骤 ...
  local MARKED_BAM="$ALIGN_DIR/${sample}.markdup.bam"; local cmd_markdup="run_picard MarkDuplicates I=\"$SORTED_BAM\" O=\"$MARKED_BAM\" M=\"$ALIGN_DIR/${sample}.markdup.metrics.txt\" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT"; run_or_skip "MARK_DUPLICATES" "$MARKED_BAM" "$cmd_markdup" "$STATUS_FILE" "$LOG_PREFIX"; local BAM_FOR_CALLING="$MARKED_BAM"
  if [[ ${#KNOWN_SITES[@]} -gt 0 ]]; then local GATK_DIR="$sample_out/gatk"; local BQSR_BAM="$GATK_DIR/${sample}.bqsr.bam"; local known_sites_args; known_sites_args=$(printf ' --known-sites %q' "${KNOWN_SITES[@]}"); local RECAL_TABLE="$GATK_DIR/${sample}.recal.table"; local cmd_bqsr="run_gatk BaseRecalibrator -R \"$REF\" -I \"$MARKED_BAM\" ${known_sites_args} -O \"$RECAL_TABLE\" && run_gatk ApplyBQSR -R \"$REF\" -I \"$MARKED_BAM\" --bqsr-recal-file \"$RECAL_TABLE\" -O \"$BQSR_BAM\" && samtools index \"$BQSR_BAM\""; run_or_skip "BQSR" "$BQSR_BAM" "$cmd_bqsr" "$STATUS_FILE" "$LOG_PREFIX"; BAM_FOR_CALLING="$BQSR_BAM"; fi
  local VARIANT_DIR="$sample_out/variants"; local GVCF="$VARIANT_DIR/${sample}.g.vcf.gz"; local RAW_VCF="$VARIANT_DIR/${sample}.raw.vcf.gz"; local FILTERED_VCF="$VARIANT_DIR/${sample}.filtered.vcf.gz"
  
  ##### 优化点 7: 在命令中动态使用 GATK HaplotypeCaller 参数 #####
  local hc_opts=""
  if [[ -n "$INTERVALS_FILE" ]]; then hc_opts="-L \"${INTERVALS_FILE}\""; fi
  local cmd_haplotypecaller="run_gatk HaplotypeCaller -R \"$REF\" -I \"$BAM_FOR_CALLING\" -O \"$GVCF\" -ERC GVCF --native-pair-hmm-threads \"$threads_to_use\" --ploidy \"${PLOIDY}\" ${hc_opts}"; 
  run_or_skip "HAPLOTYPE_CALLER" "$GVCF" "$cmd_haplotypecaller" "$STATUS_FILE" "$LOG_PREFIX"
  
  local cmd_genotypegvcfs="run_gatk GenotypeGVCFs -R \"$REF\" -V \"$GVCF\" -O \"$RAW_VCF\""; run_or_skip "GENOTYPE_GVCFS" "$RAW_VCF" "$cmd_genotypegvcfs" "$STATUS_FILE" "$LOG_PREFIX"
  
  local filter_expression="QD < ${FILTER_QD_LT} || FS > ${FILTER_FS_GT} || MQ < ${FILTER_MQ_LT} || MQRankSum < ${FILTER_MQRS_LT} || ReadPosRankSum < ${FILTER_RPRS_LT} || SOR > ${FILTER_SOR_GT} || DP < ${FILTER_DP_LT}"
  log "${LOG_PREFIX} 使用GATK过滤表达式: ${filter_expression}"
  local tmp_vcf="$VARIANT_DIR/${sample}.tmp.filtered.vcf.gz"
  local cmd_filter="run_gatk VariantFiltration -R \"$REF\" -V \"$RAW_VCF\" --filter-expression '${filter_expression}' --filter-name 'CUSTOM_FILTER' -O \"$tmp_vcf\" && bcftools view -f PASS -Oz -o \"$FILTERED_VCF\" \"$tmp_vcf\" && (command -v tabix >/dev/null 2>&1 && tabix -p vcf \"$FILTERED_VCF\") && rm -f \"${tmp_vcf}\" \"${tmp_vcf}.tbi\""
  run_or_skip "VARIANT_FILTERING" "$FILTERED_VCF" "$cmd_filter" "$STATUS_FILE" "$LOG_PREFIX"
  
  local STATS_DEPTH_FILE="$VARIANT_DIR/${sample}.ldepth.mean"; local cmd_stats="vcftools --gzvcf \"$FILTERED_VCF\" --site-mean-depth --out \"$VARIANT_DIR/${sample}\" && vcftools --gzvcf \"$FILTERED_VCF\" --site-quality --out \"$VARIANT_DIR/${sample}\""; run_or_skip "STATS_GENERATION" "$STATS_DEPTH_FILE" "$cmd_stats" "$STATUS_FILE" "$LOG_PREFIX"

  log "${LOG_PREFIX} ==== 处理完成。最终VCF文件在 $sample_out/variants 目录下 ===="
}

detect_pairs() { local dir=$1; shopt -s nullglob; local files=("$dir"/*_R1.fastq.gz "$dir"/*_R1.fq.gz "$dir"/*_1.fastq.gz "$dir"/*_1.fq.gz); shopt -u nullglob; printf '%s\n' "${files[@]}"; }
find_r2_for() { local r1=$1 r2=""; if [[ -f ${r1/_R1.fastq.gz/_R2.fastq.gz} ]]; then r2=${r1/_R1.fastq.gz/_R2.fastq.gz}; elif [[ -f ${r1/_R1.fq.gz/_R2.fq.gz} ]]; then r2=${r1/_R1.fq.gz/_R2.fq.gz}; elif [[ -f ${r1/_1.fastq.gz/_2.fastq.gz} ]]; then r2=${r1/_1.fastq.gz/_2.fastq.gz}; elif [[ -f ${r1/_1.fq.gz/_2.fq.gz} ]]; then r2=${r1/_1.fq.gz/_2.fq.gz}; fi; printf '%s' "$r2"; }
derive_sample_id() { local base; base=$(basename "$1"); local sample; sample=${base%_R1.fastq.gz}; sample=${sample%_R1.fq.gz}; sample=${sample%_1.fastq.gz}; sample=${sample%_1.fq.gz}; printf '%s' "$sample"; }

# --- 参数解析与初始化 ---
RAW_R1=""; RAW_R2=""; REF=""; OUTDIR=""; SAMPLE=""; DATA_DIR=""; KNOWN_SITES=();
PICARD_PATH="${PICARD_JAR:-}"; READ_GROUP_TEMPLATE=""; GATK_CMD="${GATK_PATH:-gatk}";
FORCE_RERUN=false; PARALLEL_JOBS=1; MIN_FREE_MEM_GB=10; CPU_FRACTION=0.8; QC_ENV_NAME=""

##### 优化点 2: 为所有新的可调参数设置默认值 #####
# fastp
FASTP_QUAL=20
FASTP_LEN=50
# BWA
BWA_USE_M=true # 默认开启-M选项
# GATK HaplotypeCaller
PLOIDY=2
INTERVALS_FILE=""
# GATK VariantFiltration
FILTER_QD_LT=1.5
FILTER_FS_GT=60.0
FILTER_MQ_LT=30.0
FILTER_MQRS_LT=-12.5
FILTER_RPRS_LT=-8.0
FILTER_SOR_GT=3.0
FILTER_DP_LT=10

##### 优化点 3: 扩展长选项解析逻辑 #####
# getopts 不支持长选项, 所以我们先手动解析它们
# 遍历所有参数，直到遇到一个不是我们定义的长选项为止
# 注意：我们把 -L 也放在这里处理，因为它更像是配置选项
extra_args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastp-qual)   FASTP_QUAL="$2"; shift 2 ;;
        --fastp-len)    FASTP_LEN="$2"; shift 2 ;;
        --bwa-use-M)    BWA_USE_M=true; shift 1 ;; # 这是一个标志，不需要值
        --ploidy)       PLOIDY="$2"; shift 2 ;;
        -L)             INTERVALS_FILE="$2"; shift 2;;
        --filter-qd)    FILTER_QD_LT="$2"; shift 2 ;;
        --filter-fs)    FILTER_FS_GT="$2"; shift 2 ;;
        --filter-mq)    FILTER_MQ_LT="$2"; shift 2 ;;
        --filter-mqrs)  FILTER_MQRS_LT="$2"; shift 2 ;;
        --filter-rprs)  FILTER_RPRS_LT="$2"; shift 2 ;;
        --filter-sor)   FILTER_SOR_GT="$2"; shift 2 ;;
        --filter-dp)    FILTER_DP_LT="$2"; shift 2 ;;
        *) # 如果不是我们定义的长选项，就收集起来留给 getopts
           extra_args+=("$1")
           shift
           ;;
    esac
done
# 将未被长选项解析器处理的参数重新放回位置参数列表
set -- "${extra_args[@]}"

##### 优化点 4: 更新 getopts 字符串，移除 -L #####
# 标准的 getopts 解析
while getopts ":1:2:r:o:s:d:k:j:m:c:q:p:g:G:fh" opt; do
  case $opt in
    1) RAW_R1=$OPTARG;; 2) RAW_R2=$OPTARG;; r) REF=$OPTARG;; o) OUTDIR=$OPTARG;;
    s) SAMPLE=$OPTARG;; d) DATA_DIR=$OPTARG;; k) KNOWN_SITES+=("$OPTARG");;
    j) PARALLEL_JOBS=$OPTARG;; m) MIN_FREE_MEM_GB=$OPTARG;;
    c) CPU_FRACTION=$OPTARG;; q) QC_ENV_NAME=$OPTARG;;
    p) PICARD_PATH=$OPTARG;; g) READ_GROUP_TEMPLATE=$OPTARG;; G) GATK_CMD=$OPTARG;;
    f) FORCE_RERUN=true;; h) usage; exit 0;;
    :) echo "错误: 参数 -$OPTARG 需要一个值" >&2; usage; exit 1;;
    \?) echo "错误: 未知参数: -$OPTARG" >&2; usage; exit 1;;
  esac
done

# ... 参数有效性与依赖检查 ...
# 首先检查conda环境
check_conda_environment

if [[ -z $REF || -z $OUTDIR ]]; then echo "错误: 缺少必填参数 -r 和 -o" >&2; usage; exit 1; fi
if [[ -n $DATA_DIR ]]; then if [[ -n $RAW_R1 || -n $RAW_R2 || -n $SAMPLE ]]; then echo "错误: 批量模式 (-d) 下不需要 -1/-2/-s" >&2; exit 1; fi; if [[ ! -d $DATA_DIR ]]; then echo "错误: 目录不存在: $DATA_DIR" >&2; exit 1; fi; else if [[ -z $RAW_R1 || -z $RAW_R2 || -z $SAMPLE ]]; then echo "错误: 单样本模式缺少必填参数 -1, -2, -s" >&2; usage; exit 1; fi; for f in "$RAW_R1" "$RAW_R2"; do if [[ ! -f $f ]]; then echo "错误: 文件不存在: $f" >&2; exit 1; fi; done; fi
if [[ ! -f $REF ]]; then echo "错误: 参考基因组不存在: $REF" >&2; exit 1; fi
##### 优化点 8: 添加对 -L 文件存在的检查 #####
if [[ -n "$INTERVALS_FILE" && ! -f "$INTERVALS_FILE" ]]; then echo "错误: 目标区域文件不存在: $INTERVALS_FILE" >&2; exit 1; fi
REQUIRED_CMDS=(fastp bwa samtools bcftools vcftools nproc bc); for cmd in "${REQUIRED_CMDS[@]}"; do if ! command -v "$cmd" >/dev/null 2>&1; then echo "错误: 缺少依赖命令: $cmd" >&2; exit 1; fi; done
if [[ -z "$QC_ENV_NAME" ]] && ! command -v fastqc >/dev/null 2>&1; then echo "错误: 缺少依赖命令: fastqc。或使用 -q 参数指定其所在的conda环境。" >&2; exit 1; fi
if [[ $GATK_CMD == gatk ]]; then if ! command -v gatk >/dev/null 2>&1; then echo "错误: 缺少依赖命令: gatk" >&2; exit 1; fi; else if [[ ! -x $GATK_CMD ]]; then echo "错误: GATK 可执行文件不可用: $GATK_CMD" >&2; exit 1; fi; fi
if [[ -n $PICARD_PATH && ! -f $PICARD_PATH ]]; then echo "错误: PICARD_JAR 文件不存在: $PICARD_PATH" >&2; exit 1; fi

# ... 线程计算 ...
TOTAL_THREADS=$(nproc); TOTAL_THREADS_TO_USE=$(printf "%.0f" "$(echo "$TOTAL_THREADS * $CPU_FRACTION" | bc)"); if (( TOTAL_THREADS_TO_USE < 1 )); then TOTAL_THREADS_TO_USE=1; fi
if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || (( PARALLEL_JOBS <= 0 )); then echo "错误: 并行任务数 (-j) 必须是一个大于0的整数" >&2; exit 1; fi
THREADS_PER_JOB=$((TOTAL_THREADS_TO_USE / PARALLEL_JOBS)); if (( THREADS_PER_JOB < 1 )); then THREADS_PER_JOB=1; fi

# ############# 主执行逻辑 #############
log "脚本启动，初始化资源配置..."; log "系统总线程数: ${TOTAL_THREADS}"; log "CPU使用率上限 (-c): ${CPU_FRACTION}"; log "流程可用总线程数: ${TOTAL_THREADS_TO_USE}"; log "并行样本数 (-j): ${PARALLEL_JOBS}"; log "将为每个样本任务分配 ${THREADS_PER_JOB} 个线程"; log "内存监控阈值 (-m): ${MIN_FREE_MEM_GB}GB"
if [[ -n "$QC_ENV_NAME" ]]; then log "FastQC 将在独立的 conda 环境中运行: $QC_ENV_NAME"; fi

log "开始全局设置：检查并创建参考基因组索引..."; if [[ ! -f "${REF}.bwt" && ! -f "${REF}.0123" ]]; then log "创建 BWA 索引..."; bwa index "$REF"; else log "BWA 索引已存在。"; fi; if [[ ! -f "${REF}.fai" ]]; then log "创建 .fai 索引..."; samtools faidx "$REF"; else log ".fai 索引已存在。"; fi; REF_DICT="${REF%.*}.dict"; if [[ ! -f "$REF_DICT" ]]; then log "创建 GATK 字典..."; run_gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"; else log "GATK 字典已存在。"; fi; log "参考基因组索引检查完毕。"

if [[ -n $DATA_DIR ]]; then
  # 批量处理模式
  mkdir -p "$OUTDIR"; mapfile -t R1_FILES < <(detect_pairs "$DATA_DIR"); if [[ ${#R1_FILES[@]} -eq 0 ]]; then echo "错误: 未在目录 $DATA_DIR 中检测到fastq文件" >&2; exit 1; fi; log "检测到 ${#R1_FILES[@]} 个样本，开始调度处理..."
  pids=(); samples_map=()
  for r1 in "${R1_FILES[@]}"; do
    while (( ${#pids[@]} >= PARALLEL_JOBS )); do
      sleep 5; finished_pids=(); running_pids=()
      for pid in "${pids[@]}"; do if ! kill -0 "$pid" 2>/dev/null; then finished_pids+=("$pid"); else running_pids+=("$pid"); fi; done
      pids=("${running_pids[@]}")
    done
    log "准备启动新任务，正在检查内存..."; check_memory_and_wait "$MIN_FREE_MEM_GB"
    r2=$(find_r2_for "$r1"); if [[ -z $r2 ]]; then log "警告: 跳过 $r1，未找到匹配的 R2 文件"; continue; fi
    sample=$(derive_sample_id "$r1"); if [[ -z $sample ]]; then log "警告: 跳过 $r1，无法解析样本 ID"; continue; fi
    sample_output_dir="$OUTDIR/$sample"; setup_dirs "$sample_output_dir"
    log "正在后台启动样本 [${sample}] 的处理... 日志: $sample_output_dir/logs/pipeline.log"
    process_sample "$sample" "$r1" "$r2" "$sample_output_dir" "$THREADS_PER_JOB" &> "$sample_output_dir/logs/pipeline.log" &
    pid=$!; pids+=("$pid"); samples_map["$pid"]="$sample"
  done
  log "所有任务均已启动，等待所有后台任务完成..."; failed_samples=(); final_exit_code=0
  for pid in "${!samples_map[@]}"; do
      sample_name=${samples_map[$pid]}
      if wait "$pid"; then log "样本 [${sample_name}] (PID: $pid) 成功完成。"; else ret_code=$?; log "!!! 错误: 样本 [${sample_name}] (PID: $pid) 处理失败，退出码: $ret_code"; failed_samples+=("$sample_name"); final_exit_code=1; fi
  done
  if (( final_exit_code == 0 )); then log "所有样本处理完毕！流程成功结束。"; else log "流程结束，但以下样本处理失败: ${failed_samples[*]}"; fi
  exit "$final_exit_code"
else
  # 单样本模式
  log "进入单样本处理模式..."; process_sample "$SAMPLE" "$RAW_R1" "$RAW_R2" "$OUTDIR" "$THREADS_PER_JOB"
fi

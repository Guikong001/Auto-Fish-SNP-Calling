#!/bin/bash

#===============================================================================
# SNP Calling Environment Setup Script
#===============================================================================
# 
# Description: Automated deployment script for SNP calling analysis environment
# Author:      wslll
# Version:     2.0.0
# Created:     $(date +%Y-%m-%d)
# License:     MIT
# 
# This script automates the setup of a complete SNP calling analysis environment
# including Conda environments, bioinformatics tools, and GATK installation.
#
# Usage:
#   ./set_env.sh [OPTIONS]
#
# Options:
#   -h, --help          Show this help message
#   -v, --verbose       Enable verbose output
#   -q, --quiet         Suppress non-essential output
#   -f, --force         Force reinstallation of existing environments
#   --skip-gatk         Skip GATK installation
#   --gatk-version      Specify GATK version (default: 4.6.2.0)
#
# Requirements:
#   - Conda/Miniconda installed and available in PATH
#   - Internet connection for downloading packages
#   - wget and unzip utilities
#
# Environments Created:
#   1. snp_call: Main analysis environment with GATK, BWA, Samtools, etc.
#   2. fastqc:   Quality control environment with FastQC
#
#===============================================================================

# Script configuration
readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_VERSION="2.0.0"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly LOG_FILE="${SCRIPT_DIR}/setup_$(date +%Y%m%d_%H%M%S).log"

# Default configuration
GATK_VERSION="4.6.2.0"
PYTHON_VERSION="3.11"
VERBOSE=false
QUIET=false
FORCE_INSTALL=false
SKIP_GATK=false

# Color codes for output formatting
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly WHITE='\033[1;37m'
readonly NC='\033[0m' # No Color

# Exit codes
readonly EXIT_SUCCESS=0
readonly EXIT_FAILURE=1
readonly EXIT_INVALID_ARGS=2

#===============================================================================
# Utility Functions
#===============================================================================

# Print formatted messages
print_header() {
    echo -e "${CYAN}===============================================================================${NC}"
    echo -e "${WHITE}$1${NC}"
    echo -e "${CYAN}===============================================================================${NC}"
}

print_step() {
    echo -e "${BLUE}>>> $1 <<<${NC}"
}

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Logging function
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
    
    case "$level" in
        "ERROR")   print_error "$message" ;;
        "WARNING") print_warning "$message" ;;
        "SUCCESS") print_success "$message" ;;
        "INFO")    [[ "$QUIET" != true ]] && print_info "$message" ;;
        *)         [[ "$QUIET" != true ]] && echo "$message" ;;
    esac
}

# Error handling
error_exit() {
    log "ERROR" "$1"
    exit "${2:-$EXIT_FAILURE}"
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Progress indicator
show_progress() {
    local duration=$1
    local message="$2"
    
    if [[ "$QUIET" == true ]]; then
        return 0
    fi
    
    echo -n "$message"
    for ((i=0; i<duration; i++)); do
        echo -n "."
        sleep 0.5
    done
    echo " 完成"
}

#===============================================================================
# Validation Functions
#===============================================================================

# Validate prerequisites
validate_prerequisites() {
    log "INFO" "验证系统先决条件..."
    
    # Check for required commands
    local required_commands=("conda" "wget" "unzip")
    for cmd in "${required_commands[@]}"; do
        if ! command_exists "$cmd"; then
            error_exit "必需的命令 '$cmd' 未找到。请确保已安装并在 PATH 中可用。"
        fi
    done
    
    # Check internet connectivity
    if ! wget --spider --quiet --timeout=10 "https://www.baidu.com" 2>/dev/null; then
        print_warning "网络连接检查失败，可能会影响包下载"
    fi
    
    # Check disk space (require at least 2GB)
    local available_space=$(df "$HOME" | awk 'NR==2 {print $4}')
    if [[ $available_space -lt 2097152 ]]; then  # 2GB in KB
        print_warning "可用磁盘空间可能不足，建议至少有 2GB 空间"
    fi
    
    log "SUCCESS" "系统先决条件验证完成"
}

# Validate GATK version format
validate_gatk_version() {
    local valid_versions=("4.6.2.0" "4.2.0.0")
    local is_valid=false
    
    for version in "${valid_versions[@]}"; do
        if [[ "$GATK_VERSION" == "$version" ]]; then
            is_valid=true
            break
        fi
    done
    
    if [[ "$is_valid" != true ]]; then
        error_exit "无效的 GATK 版本: $GATK_VERSION。支持的版本: ${valid_versions[*]}"
    fi
}

#===============================================================================
# Core Functions
#===============================================================================

# Configure Conda channels
configure_conda_channels() {
    log "INFO" "配置 Conda 软件源..."
    
    # Add channels if not already present
    local channels=("bioconda" "conda-forge")
    for channel in "${channels[@]}"; do
        if ! conda config --show channels | grep -q "$channel"; then
            conda config --add channels "$channel" || error_exit "添加 $channel 频道失败"
            log "INFO" "已添加 $channel 频道"
        else
            log "INFO" "$channel 频道已存在"
        fi
    done
    
    log "SUCCESS" "Conda 软件源配置完成"
}

# Create and setup main environment
setup_main_environment() {
    local env_name="snp_call"
    log "INFO" "设置主分析环境: $env_name"
    
    # Check if environment exists
    if conda env list | grep -q "^$env_name "; then
        if [[ "$FORCE_INSTALL" == true ]]; then
            log "WARNING" "强制重新安装，删除现有环境..."
            conda env remove -n "$env_name" -y || error_exit "删除现有环境失败"
        else
            log "WARNING" "环境 $env_name 已存在，跳过创建"
            return 0
        fi
    fi
    
    # Create environment
    log "INFO" "创建 $env_name 环境 (Python $PYTHON_VERSION)..."
    conda create -n "$env_name" python="$PYTHON_VERSION" -y || error_exit "创建环境失败"
    
    # Install MultiQC separately (known compatibility issues)
    log "INFO" "安装 MultiQC..."
    conda install -n "$env_name" -y multiqc || error_exit "MultiQC 安装失败"
    
    # Install core bioinformatics tools
    log "INFO" "安装核心生物信息学工具..."
    local packages=(
        "openjdk"
        "bwa>=0.7.17"
        "samtools>=1.17"
        "picard>=2.27"
        "fastp>=0.23"
        "bcftools>=1.17"
        "vcftools>=0.1.17"
        "seqtk>=1.3"
    )
    
    conda install -n "$env_name" -y "${packages[@]}" || error_exit "核心工具安装失败"
    
    log "SUCCESS" "主环境 $env_name 设置完成"
}

# Create FastQC environment
setup_fastqc_environment() {
    local env_name="fastqc"
    log "INFO" "设置质控环境: $env_name"
    
    # Check if environment exists
    if conda env list | grep -q "^$env_name "; then
        if [[ "$FORCE_INSTALL" == true ]]; then
            log "WARNING" "强制重新安装，删除现有环境..."
            conda env remove -n "$env_name" -y || error_exit "删除现有环境失败"
        else
            log "WARNING" "环境 $env_name 已存在，跳过创建"
            return 0
        fi
    fi
    
    # Create and setup environment
    conda create -n "$env_name" -y || error_exit "创建 FastQC 环境失败"
    conda install -n "$env_name" -c bioconda fastqc -y || error_exit "FastQC 安装失败"
    
    log "SUCCESS" "FastQC 环境设置完成"
}

# Download and install GATK
install_gatk() {
    if [[ "$SKIP_GATK" == true ]]; then
        log "INFO" "跳过 GATK 安装"
        return 0
    fi
    gatk-4.6.2.0.zip
    log "INFO" "下载并安装 GATK v$GATK_VERSION"
    
    local gatk_dir="$HOME/gatk-$GATK_VERSION"
    local gatk_zip="gatk-$GATK_VERSION.zip"
    local gatk_url="https://vip.123pan.cn/1845295995/www/Software/bioinformatics/gatk/$gatk_zip"
    
    # Change to home directory
    cd "$HOME" || error_exit "无法切换到用户主目录"
    
    # Clean up existing files if force install
    if [[ "$FORCE_INSTALL" == true ]]; then
        [[ -f "$gatk_zip" ]] && rm -f "$gatk_zip"
        [[ -d "$gatk_dir" ]] && rm -rf "$gatk_dir"
    fi
    
    # Check if GATK already exists
    if [[ -d "$gatk_dir" ]] && [[ "$FORCE_INSTALL" != true ]]; then
        log "WARNING" "GATK 目录已存在，跳过下载"
    else
        # Download GATK
        log "INFO" "从服务器下载 GATK..."
        if ! wget -O "$gatk_zip" "$gatk_url"; then
            error_exit "GATK 下载失败"
        fi
        
        # Extract GATK
        log "INFO" "解压 GATK..."
        if ! unzip -q "$gatk_zip"; then
            error_exit "GATK 解压失败"
        fi
        
        # Clean up zip file
        rm -f "$gatk_zip"
    fi
    
    log "SUCCESS" "GATK 安装完成"
}

# Configure GATK environment variables
configure_gatk_environment() {
    if [[ "$SKIP_GATK" == true ]]; then
        return 0
    fi
    
    log "INFO" "配置 GATK 环境变量"
    
    local gatk_path_cmd="export PATH=\"\$HOME/gatk-$GATK_VERSION/:\$PATH\""
    local bashrc_file="$HOME/.bashrc"
    
    # Check if already configured
    if grep -Fxq "$gatk_path_cmd" "$bashrc_file" 2>/dev/null; then
        log "INFO" "GATK 环境变量已存在于 ~/.bashrc"
    else
        log "INFO" "添加 GATK 路径到 ~/.bashrc"
        {
            echo ""
            echo "# GATK Path Configuration (Added by $SCRIPT_NAME v$SCRIPT_VERSION)"
            echo "$gatk_path_cmd"
        } >> "$bashrc_file"
        log "SUCCESS" "GATK 环境变量已写入 ~/.bashrc"
    fi
    
    # Set for current session
    export PATH="$HOME/gatk-$GATK_VERSION/:$PATH"
    log "INFO" "GATK 环境变量已在当前会话中激活"
}

# Verify installation
verify_installation() {
    log "INFO" "验证安装..."
    
    # Verify conda environments
    local environments=("snp_call" "fastqc")
    for env in "${environments[@]}"; do
        if conda env list | grep -q "^$env "; then
            log "SUCCESS" "环境 $env 验证成功"
        else
            log "ERROR" "环境 $env 验证失败"
        fi
    done
    
    # Verify GATK if not skipped
    if [[ "$SKIP_GATK" != true ]]; then
        log "INFO" "验证 GATK 安装..."
        if conda run -n snp_call gatk --version >/dev/null 2>&1; then
            local gatk_version_output=$(conda run -n snp_call gatk --version 2>&1 | head -1)
            log "SUCCESS" "GATK 验证成功: $gatk_version_output"
        else
            log "ERROR" "GATK 验证失败"
        fi
    fi
}

# Display final instructions
show_completion_message() {
    echo ""
    print_header "环境部署成功完成!"
    echo ""
    
    print_info "已创建的 Conda 环境:"
    echo "  1. snp_call: 包含 GATK 依赖(Java)、MultiQC、BWA、Samtools 等核心分析工具"
    echo "  2. fastqc:   包含 FastQC 质控工具"
    echo ""
    
    if [[ "$SKIP_GATK" != true ]]; then
        print_header "重要操作提醒"
        echo "GATK 的环境变量已配置成功。为了让它在您的终端中完全永久生效，"
        echo "请执行以下操作之一："
        echo "  1. 重新启动终端"
        echo "  2. 运行命令: source ~/.bashrc"
        echo ""
    fi
    
    print_info "使用方法:"
    echo "  激活主环境: conda activate snp_call"
    echo "  随后使用 SNP_calling 脚本开始进行自动化变异检测。"
    echo ""
    
    print_info "日志文件: $LOG_FILE"
    
    print_header "部署完成"
}

#===============================================================================
# Command Line Interface
#===============================================================================

# Show help message
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - SNP Calling Environment Setup Script

DESCRIPTION:
    Automated deployment script for SNP calling analysis environment.
    Sets up Conda environments with bioinformatics tools and GATK.

USAGE:
    $SCRIPT_NAME [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -v, --verbose           Enable verbose output
    -q, --quiet             Suppress non-essential output
    -f, --force             Force reinstallation of existing environments
    --skip-gatk             Skip GATK installation
    --gatk-version VERSION  Specify GATK version (default: $GATK_VERSION, supported: 4.6.2.0, 4.2.0.0)
    --python-version VER    Specify Python version (default: $PYTHON_VERSION)

EXAMPLES:
    $SCRIPT_NAME                    # Standard installation
    $SCRIPT_NAME --verbose          # Verbose output
    $SCRIPT_NAME --force            # Force reinstall
    $SCRIPT_NAME --skip-gatk        # Skip GATK installation
    $SCRIPT_NAME --gatk-version 4.2.0.0  # Use specific GATK version

REQUIREMENTS:
    - Conda/Miniconda installed
    - Internet connection
    - wget and unzip utilities
    - At least 2GB free disk space

For more information, visit: https://github.com/Guikong001/Auto-Fish-SNP-Calling
EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit $EXIT_SUCCESS
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -q|--quiet)
                QUIET=true
                shift
                ;;
            -f|--force)
                FORCE_INSTALL=true
                shift
                ;;
            --skip-gatk)
                SKIP_GATK=true
                shift
                ;;
            --gatk-version)
                GATK_VERSION="$2"
                validate_gatk_version
                shift 2
                ;;
            --python-version)
                PYTHON_VERSION="$2"
                shift 2
                ;;
            *)
                error_exit "未知选项: $1" $EXIT_INVALID_ARGS
                ;;
        esac
    done
}

#===============================================================================
# Main Execution
#===============================================================================

main() {
    # Initialize logging
    echo "# SNP Calling Environment Setup Log" > "$LOG_FILE"
    echo "# Started: $(date)" >> "$LOG_FILE"
    echo "# Script: $SCRIPT_NAME v$SCRIPT_VERSION" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    
    # Parse command line arguments
    parse_arguments "$@"
    
    # Set error handling
    set -eE
    trap 'error_exit "脚本在第 $LINENO 行执行失败"' ERR
    
    # Show welcome message
    print_header "欢迎使用 SNP Calling 环境自动部署脚本 v$SCRIPT_VERSION"
    
    # Validate prerequisites
    validate_prerequisites
    
    # Execute main workflow
    print_step "步骤 1/5: 配置 Conda 全局软件源"
    configure_conda_channels
    echo ""
    
    print_step "步骤 2/5: 创建并配置 Conda 环境"
    setup_main_environment
    setup_fastqc_environment
    echo ""
    
    print_step "步骤 3/5: 下载并解压 GATK"
    install_gatk
    echo ""
    
    print_step "步骤 4/5: 设置 GATK 永久环境变量"
    configure_gatk_environment
    echo ""
    
    print_step "步骤 5/5: 验证安装"
    verify_installation
    echo ""
    
    # Show completion message
    show_completion_message
    
    log "SUCCESS" "脚本执行完成"
}

# Execute main function with all arguments
main "$@"

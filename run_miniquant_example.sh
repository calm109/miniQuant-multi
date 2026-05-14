#!/bin/bash
set -e
set -o pipefail

source /home/yuyue/miniforge3/etc/profile.d/conda.sh
conda activate miniQuant-multi

cd /home/yuyue/files/miniQuant-multi
BASE_DIR="/home/yuyue/files/miniQuant-multi"
MINIQUANT_DIR="$BASE_DIR/isoform_quantification"
DATA_DIR="$BASE_DIR/example"
GTF="$DATA_DIR/annotation.gtf"
OUTPUT_DIR="$BASE_DIR/output"
mkdir -p "$OUTPUT_DIR"
THREADS=16

# ===== 日志设置开始 =====
LOG_FILE="miniQuant_run.log"
ERR_FILE="miniQuant_run.err"

# 将所有输出重定向到 LOG_FILE，将所有错误重定向到 ERR_FILE
exec > >(tee -a "$LOG_FILE") 2> >(tee -a "$ERR_FILE" >&2)
# ===== 日志设置结束 =====


# Input files from example_data
SR_SAMS=("$DATA_DIR/SR.sam" "$DATA_DIR/SR.sam")
LR_SAMS=("$DATA_DIR/LR.sam" "$DATA_DIR/LR.sam")

LR_ARGS=""
if [[ ${#LR_SAMS[@]} -gt 0 ]]; then
    LR_ARGS="-lrsam ${LR_SAMS[*]}"
fi
SR_ARGS=""
if [[ ${#SR_SAMS[@]} -gt 0 ]]; then
    SR_ARGS="-srsam ${SR_SAMS[*]}"
fi


# 1. Quantify analysis
python "$MINIQUANT_DIR/main.py" quantify \
    -gtf "$GTF" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    $LR_ARGS \
    $SR_ARGS

# 2. Identifiability analysis
python "$MINIQUANT_DIR/main.py" cal_identifiability \
    -gtf "$GTF" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    $SR_ARGS \
    $LR_ARGS

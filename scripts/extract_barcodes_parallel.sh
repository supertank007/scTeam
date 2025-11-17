#!/bin/bash
# 用法: bash grep_barcodes_parallel.sh FMD.bam NK_cells_barcodes.txt NK_cells.bam 8

bam=$1
barcode_file=$2
output_bam=$3
threads=${4:-8}

if [ $# -lt 3 ]; then
  echo "用法: $0 <input.bam> <barcodes.txt> <output.bam> [threads]"
  exit 1
fi

tmp_dir=$(mktemp -d)
tmp_sam="$tmp_dir/temp.sam"

echo "[1/4] 提取 header..."
samtools view -H "$bam" > "$tmp_sam"

echo "[2/4] 并行筛选 reads..."
cat "$barcode_file" | parallel -j "$threads" "samtools view $bam | grep {}" >> "$tmp_sam"

echo "[3/4] 转换为 BAM..."
samtools view -b "$tmp_sam" -o "$tmp_dir/tmp_unsorted.bam"

echo "[4/5] 排序并创建索引..."
samtools sort "$tmp_dir/tmp_unsorted.bam" -o "$output_bam"
samtools index "$output_bam"

echo "[5/5] 清理临时文件..."
rm -rf "$tmp_dir"

echo "✅ 完成！输出 BAM: $output_bam 以及索引: ${output_bam}.bai"

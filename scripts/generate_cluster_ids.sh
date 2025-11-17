#!/bin/bash

tmpfile=$(mktemp)
declare -A counts
celltypes=()

for file in out_cluster_*.pdf; do
  cluster=$(echo "$file" | grep -oP 'out_cluster_\K[0-9]+')
  celltype=$(echo "$file" | sed -E 's/^.*plot_//' | sed 's/.pdf$//')
  echo -e "$cluster\t$celltype"
done | sort -n | while IFS=$'\t' read -r cluster celltype; do
  counts[$celltype]=$((counts[$celltype]+1))
  count=${counts[$celltype]}
  renamed="${celltype}_${count}"
  echo "$renamed" >> "$tmpfile"
done

r_code="new_cluster_ids <- c("
r_code+=$(paste -sd "," "$tmpfile" | sed 's/\([^,]*\)/"\1"/g')
r_code+=")"

target_cfg_file="../../scoption.cfg"

tmp_cfg=$(mktemp)
grep -v '^new_cluster_ids <- c(' "$target_cfg_file" > "$tmp_cfg"
echo "$r_code" >> "$tmp_cfg"
mv "$tmp_cfg" "$target_cfg_file"

rm "$tmpfile"
echo "✅ 已成功更新 $target_cfg_file 中的 new_cluster_ids"


csv_file="../markers/top2_markers.csv"
target_cfg="../../scoption.cfg" 

tmpfile=$(mktemp)

cut -d',' -f7 "$csv_file" | tail -n +2 | sed 's/"//g' | sort -u > "$tmpfile"

r_code="markers <- c("
r_code+=$(paste -sd "," "$tmpfile" | sed 's/\([^,]*\)/"\1"/g')
r_code+=")"

tmpcfg=$(mktemp)
grep -v '^markers <- c(' "$target_cfg" > "$tmpcfg"
echo "$r_code" >> "$tmpcfg"
mv "$tmpcfg" "$target_cfg"

rm "$tmpfile"

echo "✅ 已成功更新 markers 为：$r_code"

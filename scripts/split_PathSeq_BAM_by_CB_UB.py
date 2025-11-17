# -*- coding: utf-8 -*-

import pysam
from collections import defaultdict
import os

# 创建名为 "cell" 的文件夹，如果它不存在
output_dir = "cell"
os.makedirs(output_dir, exist_ok=True)

# 打开 BAM 文件，准备读取数据
pathseq_bam = pysam.AlignmentFile("pathseq_credible_tag.bam", mode="rb")

# 创建一个默认字典来存储不同 CB 标签的读取
CB_dict = defaultdict(list)

# 遍历 BAM 文件中的每一个读取
for seg in pathseq_bam.fetch(until_eof=True):
    # 检查读取是否包含 CB 标签
    if seg.has_tag("CB"):
        CB = seg.get_tag("CB")
        CB_dict[CB].append(seg)

# 对每个 CB 标签，创建一个新的 BAM 文件，并写入相应的读取
for CB in CB_dict:
    output_path = os.path.join(output_dir, f"{CB}.bam")
    with pysam.AlignmentFile(output_path, mode="wb", template=pathseq_bam) as barcode_bam:
        for seg in CB_dict[CB]:
            barcode_bam.write(seg)

# 关闭原始 BAM 文件流，释放资源
pathseq_bam.close()

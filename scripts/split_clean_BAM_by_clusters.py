import pysam
import os

# 输入 BAM 文件
input_bam_path = "FMD.bam"
# 存放按细胞类型分割的 BAM 文件夹
output_dir = "cell_type_bam"
os.makedirs(output_dir, exist_ok=True)

# 打开原始 BAM 文件
bamfile = pysam.AlignmentFile(input_bam_path, "rb")

# 遍历每个细胞类型的 barcode 文件
barcode_files = [f for f in os.listdir(".") if f.endswith("_barcodes.txt")]

for barcode_file in barcode_files:
    cell_type = barcode_file.replace("_barcodes.txt", "")
    # 读取该细胞类型的所有 barcode
    with open(os.path.join(".", barcode_file)) as f:
        barcodes = set(line.strip() for line in f)
    
    # 输出 BAM 文件路径
    output_bam_path = os.path.join(output_dir, f"{cell_type}.bam")
    
    # 创建 BAM 文件
    with pysam.AlignmentFile(output_bam_path, "wb", template=bamfile) as out_bam:
        # 遍历原 BAM 文件
        for read in bamfile.fetch(until_eof=True):
            if read.has_tag("CB"):
                if read.get_tag("CB") in barcodes:
                    out_bam.write(read)

bamfile.close()

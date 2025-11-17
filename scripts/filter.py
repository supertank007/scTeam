# -*- coding: utf-8 -*-

import pysam
barcode_file = "list"
input_bam = "pathseq_doublefiliter_tag.bam"
output_bam = "filtered_output.bam"

with open(barcode_file, 'r') as f:
    barcodes = {line.strip() for line in f}

bamfile = pysam.AlignmentFile(input_bam, "rb")
out_bam = pysam.AlignmentFile(output_bam, "wb", template=bamfile)

for read in bamfile.fetch(until_eof=True):
    if read.query_name in barcodes:
        out_bam.write(read)
bamfile.close()
out_bam.close()


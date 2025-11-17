#!/usr/bin/env python3
import pysam
import numpy as np
import tensorflow as tf
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python script.py <input_bam>")
    sys.exit(1)

bam_input = sys.argv[1]

base_name = os.path.splitext(os.path.basename(bam_input))[0] 
output_bam = f"{base_name}_clean.bam"

model_path = "/home/code/scteamfilter.h5"
max_len = 60
batch_size = 256

print(f"Input BAM: {bam_input}")
print(f"Output BAM: {output_bam}")

model = tf.keras.models.load_model(model_path, compile=False)

bamfile = pysam.AlignmentFile(bam_input, "rb")
reads = []
names = []

for read in bamfile.fetch(until_eof=True):
    seq = read.query_sequence
    if seq is None:
        continue
    reads.append(seq[:60])      
    names.append(read.query_name)

bamfile.close()

print(f"共提取 {len(reads)} 条 reads 进行预测。")

def one_hot_encode_seq(seq, max_len=max_len):
    mapping = {'A': [1, 0, 0, 0],
               'C': [0, 1, 0, 0],
               'G': [0, 0, 1, 0],
               'T': [0, 0, 0, 1],
               'N': [0, 0, 0, 0]}
    return [mapping.get(base, [0, 0, 0, 0]) for base in seq.ljust(max_len, 'N')]

X_predict = np.array([one_hot_encode_seq(s) for s in reads])

pred_probs = model.predict(X_predict, batch_size=batch_size, verbose=1)
pred_labels = (pred_probs > 0.5).astype(int).flatten()


remove_names = set([names[i] for i, label in enumerate(pred_labels) if label == 1])
print(f"预测为污染的 reads 数量: {len(remove_names)}")

bam_in = pysam.AlignmentFile(bam_input, "rb")
bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

for read in bam_in.fetch(until_eof=True):
    if read.query_name not in remove_names:
        bam_out.write(read)

bam_in.close()
bam_out.close()
print(f"过滤完成，输出 BAM: {output_bam}")

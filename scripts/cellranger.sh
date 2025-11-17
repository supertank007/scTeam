#!/bin/bash
export PATH=/home/code/cellranger-7.0.1:$PATH
source scoption.cfg
if [ ! -f  scoption.cfg  ] ; then
        echo "scoptions.cfg not found!"
        exit

fi

if [ ! -f  scsampleList.txt  ] ; then
        echo "scsampleList.txt not found!"
        exit

fi




# ÷∂Ø…Ë÷√
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do
  if [ -f "${i}_S1_L001_I1_001.fastq.gz" ] && [ -f "${i}_S1_L001_R1_001.fastq.gz" ] && [ -f "${i}_S1_L001_R2_001.fastq.gz" ]; then
  mkdir $i
  mv ${i}_S1_L001_I1_001.fastq.gz $i/;mv ${i}_S1_L001_R1_001.fastq.gz $i/;mv ${i}_S1_L001_R2_001.fastq.gz $i/
  echo "move $i done"
  else
    if [ -f "${i}_S1_L001_R1_001.fastq.gz" ] && [ -f "${i}_S1_L001_R2_001.fastq.gz" ]; then
    mkdir $i
    mv ${i}_S1_L001_R1_001.fastq.gz $i/;mv ${i}_S1_L001_R2_001.fastq.gz $i/
    echo "move $i done"
    else 
    echo "can not found fastq files!"
    fi
  fi
done




#fastq
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do
  if [ -f "${i}_1.fastq.gz" ] && [ -f "${i}_2.fastq.gz" ] && [ -f "${i}_3.fastq.gz" ]; then
  mkdir $i
  mv ${i}_1*.gz $i/${i}_S1_L001_I1_001.fastq.gz;mv ${i}_2*.gz $i/${i}_S1_L001_R1_001.fastq.gz;mv ${i}_3*.gz $i/${i}_S1_L001_R2_001.fastq.gz
  echo "change $i name done"
  else
    if [ -f "${i}_1.fastq.gz" ] && [ -f "${i}_2.fastq.gz" ]; then
    mkdir $i
    mv ${i}_1*.gz $i/${i}_S1_L001_R1_001.fastq.gz;mv ${i}_2*.gz $i/${i}_S1_L001_R2_001.fastq.gz
    echo "change $i name done"
    else 
    echo "can not found fastq files!"
    fi
  fi
done


#fq
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do
  if [ -f "${i}_1.fq.gz" ] && [ -f "${i}_2.fq.gz" ] && [ -f "${i}_3.fq.gz" ]; then
  mkdir $i
  mv ${i}_1*.gz $i/${i}_S1_L001_I1_001.fastq.gz;mv ${i}_2*.gz $i/${i}_S1_L001_R1_001.fastq.gz;mv ${i}_3*.gz $i/${i}_S1_L001_R2_001.fastq.gz
  echo "change $i name done"
  else
    if [ -f "${i}_1.fq.gz" ] && [ -f "${i}_2.fq.gz" ]; then
    mkdir $i
    mv ${i}_1*.gz $i/${i}_S1_L001_R1_001.fastq.gz;mv ${i}_2*.gz $i/${i}_S1_L001_R2_001.fastq.gz
    echo "change $i name done"
    else 
    echo "can not found fq files!"
    fi
  fi
done




if [ $genome == "mouse" ] ; then 
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do
cellranger count --id sc_${i}  --fastqs ./$i --transcriptome /DATABANK/refdata-gex-mm10-2020-A --localcores $threads
done
fi

if [ $genome == "human" ] ; then 
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do
cellranger count --id sc_${i}  --fastqs ./$i --transcriptome /DATABANK/refdata-gex-GRCh38-2020-A --localcores $threads
done
fi

mkdir matrix
for i in `cat scsampleList.txt  |grep -v "#" |cut -f 1 `;
do 
mkdir matrix/$i
cp sc_${i}/outs/raw_feature_bc_matrix/* matrix/$i/
done
echo "You can directly analyze using the matrices in the matrix folder"

echo "all done"

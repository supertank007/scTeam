#!/bin/bash
#
#
source scoption.cfg
if [ ! -f  scoption.cfg  ] ; then
        echo "scoptions.cfg not found!"
        exit

fi

if [ ! -f  scsampleList.txt  ] ; then
        echo "scsampleList.txt not found!"
        exit

fi

  #load databank
  microbe_genome="/DATABANK/scmicrobiodata/microbio.fa"
  microbe_fasta="/DATABANK/scmicrobiodata/microbio.fa"
  microbe_fai="/DATABANK/scmicrobiodata/microbio.fa.fai"
  microbe_dict_file="/DATABANK/scmicrobiodata/microbio.dict"
  microbe_bwa_image="/DATABANK/scmicrobiodata/microbio.fa.img"
  taxonomy_db="/DATABANK/scmicrobiodata/microbio_taxonomy.db"
  contaminantfile="/DATABANK/scmicrobiodata/microbio-vecscreen-combined-matches.bed"
  gatk="/root/miniconda3/envs/gatk/bin/gatk"
##host files  
  if [ $genome == "mouse" ] ; then 
    host_hss_file="/DATABANK/refdata-gex-mm10-2020-A/fasta/genome.bfi"
    host_bwa_image="/DATABANK/refdata-gex-mm10-2020-A/fasta/genome.img"
  fi
  
  if [ $genome == "human" ] ; then 
    host_hss_file="/DATABANK/refdata-gex-GRCh38-2020-A/fasta/genome.bfi"
    host_bwa_image="/DATABANK/refdata-gex-GRCh38-2020-A/fasta/genome.img"
  fi
##
###envset
export PATH=/root/miniconda3/envs/gatk/bin:$PATH
#####

###call singlecell unaligned.bam###
#######
cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  i={};
  if [ ! -f "micro_${i}/cell_unaligned_fix.bam" ]; then
    bam_file="sc_${i}/outs/possorted_genome_bam.bam"
    mkdir -p micro_${i}
    samtools view -b -f 4 "$bam_file" > micro_${i}/unaligned.bam
    zcat matrix/${i}/barcodes.tsv.gz | sed "s/^/CB:Z:/" > micro_${i}/barcode.txt
    samtools view -h micro_${i}/unaligned.bam | grep -f micro_${i}/barcode.txt > micro_${i}/cell_unaligned.sam
    samtools view -bS micro_${i}/cell_unaligned.sam > micro_${i}/cell_unaligned.bam
    samtools view -H micro_${i}/unaligned.bam > micro_${i}/header.txt
    samtools reheader micro_${i}/header.txt micro_${i}/cell_unaligned.bam > micro_${i}/cell_unaligned_fix.bam
    rm micro_${i}/barcode.txt
    rm micro_${i}/unaligned.bam
    rm micro_${i}/cell_unaligned.sam
    rm micro_${i}/cell_unaligned.bam
    rm micro_${i}/header.txt
  else
    echo "micro_${i}/cell_unaligned_fix.bam exist"
  fi
'
###call singlecell unaligned.bam done###

###allcell pipeline start####
####
cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  i={};
    if [ ! -f "micro_${i}/cell_unaligned_final.sam" ]; then
      samtools index micro_${i}/cell_unaligned_fix.bam
      samtools fastq micro_${i}/cell_unaligned_fix.bam > micro_${i}/cell_unaligned.fastq
      trimmomatic SE -phred33 micro_${i}/cell_unaligned.fastq micro_${i}/cell_unaligned_trimmed.fastq ILLUMINACLIP:/home/code/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
      fastp -w 16 --unqualified_percent_limit 40 --cut_tail --low_complexity_filter --length_required 25 --disable_adapter_trimming -i micro_${i}/cell_unaligned_trimmed.fastq -o micro_${i}/cell_TRIMMED.fastq.gz --failed_out micro_${i}/cell_failed.fastq.gz -j micro_${i}/cell_report.json -h micro_${i}/cell_report.html
      gatk FastqToSam --FASTQ micro_${i}/cell_TRIMMED.fastq.gz --OUTPUT micro_${i}/cell_unaligned_final.sam --SM SCMICRO
      sed -i "s/!/?/g" micro_${i}/cell_unaligned_final.sam
      rm micro_${i}/cell_unaligned.fastq
      rm micro_${i}/cell_unaligned_trimmed.fastq
    else
    echo "micro_${i}/cell_unaligned_final.sam exist"
    fi
'
####

for i in `cat scsampleList.txt`;
  do
    if [ ! -f "micro_${i}/cell_pathseq.bam" ]; then
      mkdir micro_${i}/tmp
      sam_file="micro_${i}/cell_unaligned_final.sam"
      gatk --java-options "-Xmx96g -Xms96G -Djava.io.tmpdir=micro_${i}/tmp/ -XX:+UseG1GC -XX:ParallelGCThreads=16 -XX:ConcGCThreads=16" PathSeqPipelineSpark --filter-duplicates false --min-score-identity .7 --input $sam_file --filter-bwa-image $host_bwa_image --kmer-file $host_hss_file --microbe-bwa-image $microbe_bwa_image --microbe-dict $microbe_dict_file --taxonomy-file $taxonomy_db --output micro_${i}/cell_pathseq.bam --scores-output micro_${i}/cell_pathseq_output --filter-metrics filter_metrics --score-metrics score_metrics --spark-master local[32]
    else
      echo "micro_${i}/cell_pathseq.bam exist"
    fi
  done

#####Pollution removal######
for i in `cat scsampleList.txt`;
  do
    cd micro_${i}
      /root/miniconda3/envs/predict/bin/python /home/code/scteamfilter.py cell_pathseq.bam
      /root/miniconda3/envs/predict/bin/samtools sort -o sorted_cell_pathseq.bam cell_pathseq_clean.bam
      /root/miniconda3/envs/predict/bin/samtools index sorted_cell_pathseq.bam
    cd ..
  done


for i in `cat scsampleList.txt`;
  do
  bedtools intersect -abam micro_${i}/sorted_cell_pathseq.bam -b $contaminantfile > micro_${i}/pathseq.contaminants.bam
  samtools view micro_${i}/pathseq.contaminants.bam | cut -f 1 > micro_${i}/contaminants.qname.txt
  n_alignments=$(cat micro_${i}/contaminants.qname.txt | wc -l)
  if [ $n_alignments -eq 0 ]; then
      cp micro_${i}/sorted_cell_pathseq.bam micro_${i}/pathseq_credible.bam
  else
      gatk FilterSamReads \
      -I micro_${i}/sorted_cell_pathseq.bam -O micro_${i}/pathseq_credible.bam --READ_LIST_FILE micro_${i}/contaminants.qname.txt \
      --FILTER excludeReadList
  fi
  done

for i in `cat scsampleList.txt`;
  do
    gatk PathSeqScoreSpark --min-score-identity .7 --unpaired-input "micro_${i}/pathseq_credible.bam" --taxonomy-file "/DATABANK/scmicrobiodata/microbio_taxonomy.db" --scores-output "micro_${i}/credible_output"
  done

#####allcell pipeline done####
for i in `cat scsampleList.txt`;
  do
    cd "micro_${i}"
      cat cell_pathseq_output | grep "species" | awk '{print $4 "\t" $5 "\t" $8}' > tmpout
      awk -F'\t' 'BEGIN{OFS="\t"} {sum+=$3; row[NR]=$0; v[NR]=$3}
      END{for(i=1;i<=NR;i++){print row[i], v[i]/sum*100}}' tmpout > ${i}.micro
      rm tmpout
      cp ${i}.micro ../
    cd ..
  done

for i in `cat scsampleList.txt`;
  do
    cd "micro_${i}"
      cat credible_output | grep "species" | awk '{print $4 "\t" $5 "\t" $8}' > tmpout
      awk -F'\t' 'BEGIN{OFS="\t"} {sum+=$3; row[NR]=$0; v[NR]=$3}
      END{for(i=1;i<=NR;i++){print row[i], v[i]/sum*100}}' tmpout > ${i}_credible.micro
      rm tmpout
      cp ${i}_credible.micro ../
    cd ..
  done

###singlecell pipeline###
############################
###make filterlist####
ln -s /DATA/matrix/STEP3.RData /DATA/STEP3.RData
Rscript /home/code/makefilterlist.R
#####done#######
cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  cd micro_"$1"
  python /home/code/add_tags_to_PathSeq_bam.py
  python /home/code/split_PathSeq_BAM_by_CB_UB.py
' _ {}

cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  cp "barcodes_by_sample/$1.barcodes.txt" "micro_$1/list"
    cd micro_"$1"
    mkdir cellfilter
      for i in `cat list`;do mv cell/${i}.bam cellfilter/ ;done
    rm -r cell/
    mkdir Pathseq
    ls cellfilter > tmplist
      for n in `cat tmplist`;do gatk PathSeqScoreSpark --min-score-identity .5 --unpaired-input "cellfilter/$n" --taxonomy-file "/DATABANK/scmicrobiodata/microbio_taxonomy.db" --scores-output "Pathseq/${n}.txt" ;done
    echo "PATHSEQ COMPLETE"
' _ {}

######PATHSEQ COMPLETE######
cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  sample="$1"
  cd "micro_${sample}" || exit 1
  input_dir="Pathseq"
  output_dir="cellR"
  mkdir -p "$output_dir"

  for file in "$input_dir"/*; do
    if [ -f "$file" ]; then
      line_count=$(wc -l < "$file")
      if [ "$line_count" -gt 2 ]; then
        cp "$file" "$output_dir/"
      fi
    fi
  done
' _ {}
######merge socore######







######START allcell KRAKEN2##########
cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  sample="$1"
  cd "micro_${sample}" || exit 1
  echo "START KRAKEN2"
  mkdir KRAKEN2
  mkdir KRAKEN2/allcell
    kraken2 --db "/DATABANK/minikraken_8GB_20200312/" "cell_TRIMMED.fastq.gz" --output KRAKEN2/allcell/${sample}_output.txt --report KRAKEN2/allcell/${sample}_report.txt
    cat KRAKEN2/allcell/${sample}_report.txt |grep S1|awk -F "\t" '{print $6"\t"$1"\t"$2}'|sed 's/^[ \t]*\([A-Za-z]\)/\1/'|sed 's/ /_/g'|awk 'BEGIN {sum=0} {sum += $3; data[NR] = $0} END {for (i=1; i<=NR; i++) {split(data[i], a); print a[1]"\t" a[2]"\t"a[3]"\t" a[3]/sum}}' > KRAKEN2/allcell/${sample}.speic
' _ {}  



cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  sample="$1"
  cd "micro_${sample}"
  echo "START KRAKEN2"

  mkdir -p KRAKEN2/allcell

  kraken2 --db "/DATABANK/minikraken_8GB_20200312/" \
          "cell_TRIMMED.fastq.gz" \
          --output KRAKEN2/allcell/${sample}_output.txt \
          --report KRAKEN2/allcell/${sample}_report.txt

  cat KRAKEN2/allcell/${sample}_report.txt | grep S1 | \
    awk -F "\t" '"'"'{print $6"\t"$1"\t"$2}'"'"' | \
    sed '"'"'s/^[ \t]*\([A-Za-z]\)/\1/'"'"' | \
    sed '"'"'s/ /_/g'"'"' | \
    awk '"'"'BEGIN {sum=0} {sum += $3; data[NR] = $0} \
          END {for (i=1; i<=NR; i++) {split(data[i], a); \
          print a[1]"\t" a[2]"\t"a[3]"\t" a[3]/sum}}'"'"' \
    > KRAKEN2/allcell/${sample}.speic
' _ {}


############################

cat scsampleList.txt | xargs -P 8 -I {} bash -c '
  sample="$1"
  cd "micro_${sample}" || exit 1
  echo "START KRAKEN2"
  mkdir KRAKEN2
  mkdir KRAKEN2/cell
  mkdir KRAKEN2/cellfilter
  mkdir KRAKEN2/cellfastq
  python /home/code/s
  for i in `cat list`;do mv KRAKEN2/cell/${i}.bam KRAKEN2/cellfilter/ done
  rm -r mkdir KRAKEN2/cell
  ls KRAKEN2/cellfilter > KRAKEN2/tmplist
  for i in `cat KRAKEN2/tmplist|sed 's/.bam//g'`;do samtools fastq KRAKEN2/cellfilter/${i}.bam > KRAKEN2/cellfastq/${i}.fastq done
  for i in `cat KRAKEN2/tmplist|sed 's/.bam//g'`;do kraken2 --db "/DATABANK/minikraken_8GB_20200312/" "KRAKEN2/cellfastq/${i}.fastq" --output KRAKEN2/${i}_output.txt --report KRAKEN2/${i}_report.txt done
  
  

gatk PathSeqScoreSpark --min-score-identity .7 --unpaired-input "GGCCGATGTCCGCTGA-1.bam" --taxonomy-file "/DATABANK/scmicrobiodata/microbio_taxonomy.db" --scores-output "sososo.txt"



awk 'BEGIN {sum=0} {sum += $3; data[NR] = $0} END {for (i=1; i<=NR; i++) {split(data[i], a); print a[1], a[2], a[3], a[3]/sum}}'

 cat report.txt |grep S1|awk -F "\t" '{print $6"\t"$1"\t"$2}'|awk 'BEGIN {sum=0} {sum += $3; data[NR] = $0} END {for (i=1; i<=NR; i++) {split(data[i], a); print a[1]"\t" a[2]"\t"a[3]"\t" a[3]/sum}}'


cat SRR28551839_report.txt |grep S1|awk -F "\t" '{print $6"\t"$1"\t"$2}'|sed 's/^[ \t]*\([A-Za-z]\)/\1/'|sed 's/ /_/g'|awk 'BEGIN {sum=0} {sum += $3; data[NR] = $0} END {for (i=1; i<=NR; i++) {split(data[i], a); print a[1]"\t" a[2]"\t"a[3]"\t" a[3]/sum}}' > SRR28551839.speic
#!/bin/bash
if [ -d "/DATA/matrix" ]; then
    cd /DATA/matrix
    echo "PWD set /DATA/matrix"
else
    cd /DATA
    echo "PWD set /DATA"
fi

source scoption.cfg
  if [ ! -f  scoption.cfg  ] ; then
	  echo "scoptions.cfg not found!"
	  exit
  fi

###STEP1###

if [ "$automatic_resolution" == "1" ]; then
    echo "use scRNAstep1_auto_resolution.R"
    Rscript /home/code/scRNAstep1_auto_resolution.R
elif [ "$manual_resolution" == "1" ]; then
    echo "use scRNAstep1_manual_resolution.R"
    Rscript /home/code/scRNAstep1_manual_resolution.R
else
    echo "please set automatic_resolution or manual_resolution = 1"
fi

###STEP2###
if [ "$automatic_resolution" == "1" ]; then
    Rscript /home/code/scRNAstep2_clusters.R
elif [[ "$manual_resolution" == "1" && -n "$resolution" ]]; then
    echo "? manual_resolution = $resolution"
    Rscript /home/code/scRNAstep2_clusters.R
    Rscript /home/code/scRNAstep3_SingleR.R
else
    echo "No manual parameters are set, automatic process is used"
fi
###STEP3###

if [ "$automatic_annotation" == "1" ]; then
  mkdir results/AIANNO
    mv cluster_*.csv results/AIANNO
    cd results/AIANNO
    for i in `ls cluster_*.csv`;
      do python /home/code/celltype.py $i Adult_Mouse_Gut.pkl
         python /home/code/plot.py out_${i}.csv
      done
    
    
elif[ "$manual_annotation" == "1" ]; then
      Rscript /home/code/scRNAstep4_annotation.R
else
  echo "No manual parameters are set, automatic process is used"
fi

########scMicro#######

######CPDB######
      if [ $cpdb -eq 1 ]; then
  Rscript /home/code/scRNAstep5_cellphonedb.R
  cellphonedb method statistical_analysis cpdb_meta.txt cpdb_counts.txt --threads $threads
  cellphonedb plot dot_plot
  cellphonedb plot heatmap_plot cpdb_meta.txt
  Rscript /home/code/scRNAstep5_cellphonedb_plot.R
      else
  echo "Not conducting cell communication analysis,Thank you good luck"
      fi
###
	if [ $KEGGGO -eq 1 ]; then
	echo "Start Enrich Analysis"
	 Rscript /home/code/scRNAstep7_KEGGGO.R
	echo "KEGG GO enrich complete!"
	fi
	
	if [ $monocle -eq 1 ]; then
        echo "Start monocle Analysis"
         Rscript /home/code/scRNAstep8_monocle.R
        echo "monocle Analysis complete!"
	fi
###
      if [ $TF -eq 1 ]; then
  echo "Start Transcription Factor Analysis"
  mkdir scenic
  python /usr/bin/scTFchange.py
            if [ "$genome" = "mouse" ]; then
    echo "mouse genome"
    #grn
    pyscenic grn \
      --num_workers $threads \
      --output ./scenic/adj.sample.tsv \
      --method grnboost2 \
      ./scenic/sample.loom \
      /DATABANK/mm_mgi_tfs.txt
    #ctx
    pyscenic ctx \
    ./scenic/adj.sample.tsv \
    /DATABANK/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname /DATABANK/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname ./scenic/sample.loom \
    --mode "dask_multiprocessing" \
    --output ./scenic/reg.csv \
    --num_workers $threads \
    --mask_dropouts
    #aucell
    pyscenic aucell \
    ./scenic/sample.loom \
    ./scenic/reg.csv \
    --output ./scenic/sample_SCENIC.loom \
    --num_workers $threads
    
    Rscript /home/code/scRNAstep6_scenic.R
            else
    echo "human genome"
    #grn
    pyscenic grn \
      --num_workers $threads \
      --output ./scenic/adj.sample.tsv \
      --method grnboost2 \
      ./scenic/sample.loom \
      /DATABANK/hs_hgnc_tfs.txt
    #ctx
    pyscenic ctx \
    ./scenic/adj.sample.tsv \
    /DATABANK/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname /DATABANK/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname ./scenic/sample.loom \
    --mode "dask_multiprocessing" \
    --output ./scenic/reg.csv \
    --num_workers $threads \
    --mask_dropouts
    #aucell
    pyscenic aucell \
    ./scenic/sample.loom \
    ./scenic/reg.csv \
    --output ./scenic/sample_SCENIC.loom \
    --num_workers $threads
    
    Rscript /home/code/scRNAstep6_scenic.R
            fi
      else
      echo "Do not perform transcription factor analysis,Thank you good luck"
      fi
        else echo  "without manual annotation"
        fi
	fi
        echo "all done! good luck!"
fi

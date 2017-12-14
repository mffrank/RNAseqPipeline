#!/bin/bash

# Iterate over the directories with the raw files and run the pipeline
# Goto folder containing the subfolders with fastq files and run script
for i in *;
do
  cd $i;
  /nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/HeLa_Celllines/results/RNAseq_pipeline/rnaseq_pipeline.sh 0 0 $i /nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/HeLa_Celllines/data/RNAseq/Alignment/$i /nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/HeLa_Celllines/results/RNAseq_pipeline/rnaseq_pipeline_config.sh *.fastq;
  cd ..;
done

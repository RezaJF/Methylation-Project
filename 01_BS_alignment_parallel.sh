#!/bin/bash
#SBATCH -p unlimited
#SBATCH --job-name=01_BS_alignment_pipeline
#SBATCH -o 01_BS_OFL_B_alignment_pipeline.out
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32GB
#SBATCH -t 120:00:00

source ~/.bashrc

#conda activate biobase

#######################
#include variables from file called "variables"
#######################

source variables

#source ~/.bs_mapping.sh

export LC_ALL=C
export MALLOC_ARENA_MAX=4 #Sets the maximum number of memory pools used, regardless of the number of cores

#mkdir ../processing
#mkdir ../processing/1-fastqc
#mkdir ../processing/2-trim
#mkdir ../processing/3-fastqc
#mkdir ../processing/4-bismark

##########################
#Concatentate Multilane FQ Files if required
##########################
lane1_pair1=$1
lane2_pair1=$2
lane3_pair1=$3

lane1_pair2=$4
lane2_pair2=$5
lane3_pair2=$6

sampleID=($(echo $lane1_pair1 | cut -d "_" -f 1))

cat ${FASTQ_PATH}/$lane1_pair1 ${FASTQ_PATH}/$lane2_pair1 ${FASTQ_PATH}/$lane3_pair1 > "$sampleID"_cat1.fq.gz
cat ${FASTQ_PATH}/$lane1_pair2 ${FASTQ_PATH}/$lane2_pair2 ${FASTQ_PATH}/$lane3_pair2 > "$sampleID"_cat2.fq.gz


#cat FASTQ_list | while read lane1_pair1 lane1_pair2;
#  do
#    if [ $no_lanes -eq "1" ]; then
      fq_align_1=$lane1_pair1
      fq_align_2=$lane1_pair2

#    fi

#      if [ $no_lanes -eq "2"  ]; then
#        cat $lane1_pair1 $lane2_pair1 > "$sampleID"_cat1.fq
#        cat $lane1_pair2 $lane2_pair2 > "$sampleID"_cat2.fq
#        fq_align_1="$sampleID"_cat1.fq
#        fq_align_2="$sampleID"_cat2.fq
#    fi

#      if [ $no_lanes -eq "3" ]; then
#        cat $lane1_pair1 $lane2_pair1 $lane3_pair1 > "$sampleID"_cat1.fq
#        cat $lane1_pair2 $lane2_pair2 $lane3_pair3 > "$sampleID"_cat2.fq
#        fq_align_1="$sampleID"_cat1.fq
#        fq_align_2="$sampleID"_cat2.fq

#    fi


#      sampleID=($(echo $lane1_pair1 | cut -d "/" -f 8 | cut -d "_" -f 1,2))

 
	echo $sampleID



      # 1. fastqc
      echo "fastqc 1st start on"; date; echo "====="


      fastqc -o $PATH_DIR/2nd_Run/processing/1-fastqc -f fastq "$sampleID"_cat1.fq.gz "$sampleID"_cat2.fq.gz

      echo "====="; echo "Pre-trim fastqc analysis completed on"; date

      # 2. trim galore, manual correct rrbs enzyme parameters
      echo "trim galore start on"; date; echo "====="



      trim_galore -o $PATH_DIR/2nd_Run/processing/2-trim \
                  --paired \
		  --gzip \
                  --retain_unpaired \
		  --clip_R1 15 \
		  --clip_R2 15 \
                  --three_prime_clip_R1 15 \
                  --three_prime_clip_R2 15 \
                  --basename $sampleID \
                  "$sampleID"_cat1.fq.gz \
                  "$sampleID"_cat2.fq.gz

      # 3. fastqc
      echo "Post-trim fastqc start on"; date; echo "====="


      fastqc -o $PATH_DIR/2nd_Run/processing/3-fastqc \
              -f fastq $PATH_DIR/2nd_Run/processing/2-trim/${sampleID}_val_1.fq.gz \
              $PATH_DIR/2nd_Run/processing/2-trim/${sampleID}_val_2.fq.gz

      echo "====="; echo "Post-trim fastqc end on"; date

      # 4.1 bismark
      echo "bismark start on"; date; echo "====="

      bismark --bowtie2  -N 1 -p 2 --non_directional --bam \
      -o $PATH_DIR/2nd_Run/processing/4-bismark \
      --temp_dir $PATH_DIR/2nd_Run/processing/4-bismark \
      $refgenome \
      --parallel 4 \
      -1 $PATH_DIR/2nd_Run/processing/2-trim/${sampleID}_val_1.fq.gz \
      -2 $PATH_DIR/2nd_Run/processing/2-trim/${sampleID}_val_2.fq.gz


      # 4.2  deduplicate_bismark
      echo "deduplicate_bismark start on"; date; echo "====="

      deduplicate_bismark -p --bam $PATH_DIR/2nd_Run/processing/4-bismark/${sampleID}_val*_bismark_bt2_pe.bam

      echo "====="; echo "deduplicate_bismark end on"; date

      # 4.3 bismark extractor
      echo "bismark methylation extractor start on"; date; echo "====="

      bismark_methylation_extractor --comprehensive \
                                    --merge_non_CpG \
                                    --report \
                                    --paired-end \
                                    -p --no_overlap \
                                    ${sampleID}_val*_bismark_bt2_pe.deduplicated.bam \
                                    -o $PATH_DIR/2nd_Run/processing/6-methylation_extraction

      echo "====="; echo "bismark methylation extractor end on"; date


#done

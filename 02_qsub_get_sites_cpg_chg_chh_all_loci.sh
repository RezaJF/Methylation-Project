#!/bin/bash

conda activate python2

script_PATH=/Users/reza/Documents/methylation_project/script_bismark

cat sample_list.txt | while read sample;
  do

    samplename=($(echo $sample | cut -d "_" -f 3))

    # CpG Context:

    awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' $sample | sort -k1,1d -k2,2n > ${samplename}_cpg_sorted.txt
    python ${script_PATH}/bs_MergePeReads.py -i ${samplename}_cpg_sorted.txt -o ${samplename}_cpg.txt

    rm ${samplename}_cpg_sorted.txt


    echo $samplename "processed at:";date

  done

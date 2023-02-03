#!/bin/bash

conda activate python2

script_PATH=/Users/reza/Documents/methylation_project/script_bismark
refconfig=/Users/reza/Documents/methylation_project/script_bismark/hs_GRCh37.73_config.txt
file_PATH=/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/02_CpG_context

### Calculate 3k window methylation frequency for each sample

cat SC_bulk_list.txt | while read SC bulk SC_baseName Bulk_baseNeme;
  do
    grep -v "Bismark" ${file_PATH}/$SC | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > ${SC_baseName}_cpg_sorted.txt
    python ${script_PATH}/bs_MergePeReads.py -i ${SC_baseName}_cpg_sorted.txt -o ${SC_baseName}_cpg_sorted_numerical.txt
    python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i ${SC_baseName}_cpg_sorted_numerical.txt -o ${SC_baseName}_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

    sed '1d' ${SC_baseName}_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > ${SC_baseName}_cpg_bin_3k_cpg_bin_3k_sorted.txt

    rm ${SC_baseName}_cpg_sorted.txt ${SC_baseName}_cpg_sorted_numerical.txt ${SC_baseName}_cpg_bin_3k.txt

  done

  cat SC_bulk_list.txt | while read SC bulk SC_baseName Bulk_baseNeme;

    do

      paste ${Bulk_baseNeme}_cpg_bin_3k_cpg_bin_3k_sorted.txt ${SC_baseName}_cpg_bin_3k_cpg_bin_3k_sorted.txt > ${Bulk_baseNeme}-${SC_baseName}_cpg_bin_3k_sorted.txt
      awk 'length($1)<=2' ${Bulk_baseNeme}-${SC_baseName}_cpg_bin_3k_sorted.txt | awk '$7>=5 && $16>=5' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $14 "\t" $15}'  > ${Bulk_baseNeme}-${SC_baseName}_cpg_bin_3k_5cpgs_sorted.txt

      Rscript ${script_PATH}/calculate_pooled_standard_error.R ${Bulk_baseNeme} ${SC_baseName} ${Bulk_baseNeme}-${SC_baseName}.dmr.txt
      Rscript ${script_PATH}/RJF_calculate_variance_noise_mouse.R ${Bulk_baseNeme} ${SC_baseName} ${Bulk_baseNeme}-${SC_baseName}_variance_noise.txt

done

mkdir DMR && mv dmr* DMR

for i in $(ls *variance_noise.txt);
  do
    echo $i | cut -d "_" -f 1,2,3 >> x
    cat $i >> y
done
paste x y > IR_TCells_noise.txt
rm x y

#!/bin/bash

conda activate python2

#######################################################################################
####### THIS SCRIPT CALCULATE CORRELATION BETWEEN BULK & PSEUDO-BULK SAMPLES ##########
#######################################################################################


script_PATH=/Users/reza/Documents/methylation_project/script_bismark
refconfig=/Users/reza/Documents/methylation_project/script_bismark/hs_GRCh37.73_config.txt

###########################################################
### 1) concatanating single-cells to pesuedo bulk  ########
###########################################################

## 1.1) CFC batch ##
cat $(cat CFC_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_CFC_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_CFC_single_cells.txt | sort -k1,1d -k2,2n > merged_CFC_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_CFC_pseudo_bulk_cpg_sorted.txt -o CFC_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i CFC_pseudo_bulk_cpg.txt -o merged_CFC_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_CFC_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_CFC_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_CFC_single_cells.txt merged_CFC_pseudo_bulk_cpg_sorted.txt CFC_pseudo_bulk_cpg.txt merged_CFC_pseudo_Bulk_cpg_bin_3k.txt

mv merged_CFC_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk


## 1.2) CFL batch ##
cat $(cat CFL_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_CFL_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_CFL_single_cells.txt | sort -k1,1d -k2,2n > merged_CFL_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_CFL_pseudo_bulk_cpg_sorted.txt -o CFL_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i CFL_pseudo_bulk_cpg.txt -o merged_CFL_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_CFL_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_CFL_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_CFL_single_cells.txt merged_CFL_pseudo_bulk_cpg_sorted.txt CFL_pseudo_bulk_cpg.txt merged_CFL_pseudo_Bulk_cpg_bin_3k.txt

mv merged_CFL_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 1.3) OFC batch ##
cat $(cat OFC_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_OFC_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_OFC_single_cells.txt | sort -k1,1d -k2,2n > merged_OFC_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_OFC_pseudo_bulk_cpg_sorted.txt -o OFC_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i OFC_pseudo_bulk_cpg.txt -o merged_OFC_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_OFC_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_OFC_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_OFC_single_cells.txt merged_OFC_pseudo_bulk_cpg_sorted.txt OFC_pseudo_bulk_cpg.txt merged_OFC_pseudo_Bulk_cpg_bin_3k.txt

mv merged_OFC_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 1.4) OFL batch ##
cat $(cat OFL_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_OFL_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_OFL_single_cells.txt | sort -k1,1d -k2,2n > merged_OFL_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_OFL_pseudo_bulk_cpg_sorted.txt -o OFL_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i OFL_pseudo_bulk_cpg.txt -o merged_OFL_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_OFL_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_OFL_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_OFL_single_cells.txt merged_OFL_pseudo_bulk_cpg_sorted.txt OFL_pseudo_bulk_cpg.txt merged_OFL_pseudo_Bulk_cpg_bin_3k.txt

mv merged_OFL_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 1.5) YFC batch ##
cat $(cat YFC_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_YFC_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_YFC_single_cells.txt | sort -k1,1d -k2,2n > merged_YFC_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_YFC_pseudo_bulk_cpg_sorted.txt -o YFC_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i YFC_pseudo_bulk_cpg.txt -o merged_YFC_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_YFC_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_YFC_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_YFC_single_cells.txt merged_YFC_pseudo_bulk_cpg_sorted.txt YFC_pseudo_bulk_cpg.txt merged_YFC_pseudo_Bulk_cpg_bin_3k.txt

mv merged_YFC_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 1.6) YFL batch ##
cat $(cat YFL_single_cells.txt) | grep -v "Bismark" > CpG_context_merged_YFL_single_cells.txt
awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' CpG_context_merged_YFL_single_cells.txt | sort -k1,1d -k2,2n > merged_YFL_pseudo_bulk_cpg_sorted.txt

python ${script_PATH}/bs_MergePeReads.py -i merged_YFL_pseudo_bulk_cpg_sorted.txt -o YFL_pseudo_bulk_cpg.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i YFL_pseudo_bulk_cpg.txt -o merged_YFL_pseudo_Bulk_cpg_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' merged_YFL_pseudo_Bulk_cpg_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_YFL_pseudo_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_merged_YFL_single_cells.txt merged_YFL_pseudo_bulk_cpg_sorted.txt YFL_pseudo_bulk_cpg.txt merged_YFL_pseudo_Bulk_cpg_bin_3k.txt

mv merged_YFL_pseudo_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

#########################################################################
########################## 2) Parsing bulk data #########################
#########################################################################

## 2.1) CFC batch ##
sed '1d' CpG_context_CFCB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_CFC_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_CFC_Bulk_sorted.txt -o CFC_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i CFC_Bulk_CpG.txt -o CFC_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' CFC_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > CFC_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_CFC_Bulk_sorted.txt CFC_Bulk_CpG.txt CFC_bulk_CpG_bin_3k.txt
mv CFC_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk


## 2.2) CFL batch ##
sed '1d' CpG_context_CFLB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_CFL_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_CFL_Bulk_sorted.txt -o CFL_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i CFL_Bulk_CpG.txt -o CFL_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' CFL_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > CFL_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_CFL_Bulk_sorted.txt CFL_Bulk_CpG.txt CFL_bulk_CpG_bin_3k.txt
mv CFL_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk


## 2.3) OFC batch ##
sed '1d' CpG_context_OFCB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_OFC_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_OFC_Bulk_sorted.txt -o OFC_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i OFC_Bulk_CpG.txt -o OFC_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' OFC_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > OFC_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_OFC_Bulk_sorted.txt OFC_Bulk_CpG.txt OFC_bulk_CpG_bin_3k.txt
mv OFC_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 2.4) OFL batch ##
sed '1d' CpG_context_OFLB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_OFL_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_OFL_Bulk_sorted.txt -o OFL_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i OFL_Bulk_CpG.txt -o OFL_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' OFL_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > OFL_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_OFL_Bulk_sorted.txt OFL_Bulk_CpG.txt OFL_bulk_CpG_bin_3k.txt
mv OFL_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 2.5) YFC batch ##
sed '1d' CpG_context_YFCB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_YFC_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_YFC_Bulk_sorted.txt -o YFC_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i YFC_Bulk_CpG.txt -o YFC_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' YFC_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > YFC_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_YFC_Bulk_sorted.txt YFC_Bulk_CpG.txt YFC_bulk_CpG_bin_3k.txt
mv YFC_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk


## 2.5) YFL batch ##
sed '1d' CpG_context_YFLB_val_1_bismark_bt2_pe.deduplicated.txt | awk '{print $3 "\t" $4 "\t" $2 "\t" $5}' | sort -k1,1d -k2,2n > CpG_context_YFL_Bulk_sorted.txt
python ${script_PATH}/bs_MergePeReads.py -i CpG_context_YFL_Bulk_sorted.txt -o YFL_Bulk_CpG.txt
python ${script_PATH}/bs_window_merge_whole_genome_v2.0.py -i YFL_Bulk_CpG.txt -o YFL_bulk_CpG_bin_3k.txt -w 3000 -s 600 -f $refconfig

sed '1d' YFL_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > YFL_Bulk_CpG_bin_3k_sorted.txt

rm CpG_context_YFL_Bulk_sorted.txt YFL_Bulk_CpG.txt YFL_bulk_CpG_bin_3k.txt
mv YFL_Bulk_CpG_bin_3k_sorted.txt /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

#############################################################################
##################### 3) Merge buld & pseudo-bulk  ##########################
#############################################################################
cd /Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk

## 3.1) CFC batch ##
paste merged_CFC_pseudo_Bulk_CpG_bin_3k_sorted.txt CFC_Bulk_CpG_bin_3k_sorted.txt > merged_CFC_pseudoBulk_and_bulk_CpG_bin_3k.txt

## 3.2) CFL batch ##
paste merged_CFL_pseudo_Bulk_CpG_bin_3k_sorted.txt CFL_Bulk_CpG_bin_3k_sorted.txt > merged_CFL_pseudoBulk_and_bulk_CpG_bin_3k.txt

## 3.3) OFC batch ##
paste merged_OFC_pseudo_Bulk_CpG_bin_3k_sorted.txt OFC_Bulk_CpG_bin_3k_sorted.txt > merged_OFC_pseudoBulk_and_bulk_CpG_bin_3k.txt

## 3.4) OFL batch ##
paste merged_OFL_pseudo_Bulk_CpG_bin_3k_sorted.txt OFL_Bulk_CpG_bin_3k_sorted.txt > merged_OFL_pseudoBulk_and_bulk_CpG_bin_3k.txt

## 3.5) YFC batch ##
paste merged_YFC_pseudo_Bulk_CpG_bin_3k_sorted.txt YFC_Bulk_CpG_bin_3k_sorted.txt > merged_YFC_pseudoBulk_and_bulk_CpG_bin_3k.txt

## 3.5) YFL batch ##
paste merged_YFL_pseudo_Bulk_CpG_bin_3k_sorted.txt YFL_Bulk_CpG_bin_3k_sorted.txt > merged_YFL_pseudoBulk_and_bulk_CpG_bin_3k.txt

####################################################################
#################### 4) Methylation frequency  #####################
####################################################################

mkdir tmp

# 4.1 Remove aletenative assembely contigs:
# CFC
awk 'length($1)<=2' merged_CFC_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_CFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_CFC_pseudoBulk_and_bulk_CpG_bin_3k.text
# CFL
awk 'length($1)<=2' merged_CFL_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_CFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_CFL_pseudoBulk_and_bulk_CpG_bin_3k.text
# OFC
awk 'length($1)<=2' merged_OFC_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_OFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_OFC_pseudoBulk_and_bulk_CpG_bin_3k.text
# OFL
awk 'length($1)<=2' merged_OFL_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_OFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_OFL_pseudoBulk_and_bulk_CpG_bin_3k.text
# YFC
awk 'length($1)<=2' merged_YFC_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_YFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_YFC_pseudoBulk_and_bulk_CpG_bin_3k.text
# YFL
awk 'length($1)<=2' merged_YFL_pseudoBulk_and_bulk_CpG_bin_3k.txt | sort -k1,1 -k2,2n -T ./tmp | uniq > merged_YFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt
rm merged_YFL_pseudoBulk_and_bulk_CpG_bin_3k.text



# 4.2 Filter the cpgsites less than/equal 5 (retain windows with cpg >=5):
awk '$7>=5 && $16>=5' merged_CFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_CFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt
awk '$7>=5 && $16>=5' merged_CFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_CFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt
awk '$7>=5 && $16>=5' merged_OFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_OFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt
awk '$7>=5 && $16>=5' merged_OFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_OFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt
awk '$7>=5 && $16>=5' merged_YFC_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_YFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt
awk '$7>=5 && $16>=5' merged_YFL_pseudoBulk_and_bulk_CpG_bin_3k_Homo_sapiens_sorted.txt | awk '{print $9 "\t" $18}' > merged_YFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt

#!/bin/bash
#SBATCH -p unlimited
#SBATCH --job-name=MethylExtractor_auxiliary
#SBATCH -o 01.1_MethylExtractor_auxiliary.out
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16GB
#SBATCH -t 48:10:00

#source ~/.bashrc
source variables

cat FASTQ_list | while read lane1_pair1 lane1_pair2;
  do
	  sampleID=($(echo $lane1_pair1 | cut -d "_" -f 1))
	  echo $sampleID

	  echo "bismark methylation extractor for sample" $sampleID "start on"; date; echo "====="
	 
	  bismark_methylation_extractor --comprehensive \
		  		    --merge_non_CpG \
                                    --bedGraph \
                                    --buffer_size 10G \
                                    --cytosine_report \
                                    --report \
                                    --paired-end \
				    --genome_folder $refgenome \
                                    -p --no_overlap \
				    /gs/gsfs0/users/rjabalamel/methylation_project/IR_TCells/2nd_Run/processing/5-deduplication/${sampleID}_val*_bismark_bt2_pe.deduplicated.bam \
				    -o /gs/gsfs0/users/rjabalamel/methylation_project/IR_TCells/2nd_Run/processing/7-GW_cytosine_report

      echo "====="; echo "bismark methylation extractor end on"; date


done

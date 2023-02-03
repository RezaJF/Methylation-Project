#!/bin/bash

cat FASTQ_list | while read lane1_pair1 lane2_pair1 lane3_pair1 lane1_pair2 lane2_pair2 lane3_pair2;

	do
		sbatch 01_BS_alignment_parallel.sh $lane1_pair1 $lane2_pair1 $lane3_pair1 $lane1_pair2 $lane2_pair2 $lane3_pair2

	done	

#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N GB1_D2regA2S1T1
#$ -l h_rt=6:00:00,h_vmem=8G,virtual_free=8G
#$ -t 1-10
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_GB1_dG_estimation.R \
		--dataset_name GB1_dG_dataset2reg \
		-a 2 \
		-i ${SGE_TASK_ID}
	

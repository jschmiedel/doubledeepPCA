#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N GB1_4sD2A2S1T1
#$ -l h_rt=6:00:00,h_vmem=8G,virtual_free=8G
#$ -t 11-100
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_GB1_dG_4state_estimation.R \
		--dataset_name GB1_dG_dataset2 \
		-a 2 \
		-i ${SGE_TASK_ID}
	

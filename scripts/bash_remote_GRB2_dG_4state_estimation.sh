#!/bin/bash
#$ -q long-sl7
#$ -N GRB2_4s_2_1e9
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l h_rt=100:00:00,h_vmem=55G,virtual_free=55G
#$ -pe smp 15
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_GRB2_dG_4state_estimation.R \
		-m 2 \
		-b 0
	

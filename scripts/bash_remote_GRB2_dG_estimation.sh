#!/bin/bash
#$ -q long-sl7
#$ -N dG_1_bf0
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l h_rt=100:00:00
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -t 1-10

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_GRB2_dG_estimation.R \
	-s ${SGE_TASK_ID} \
	-m 1 \
	-d GRB2_dG_dataset

#!/bin/bash
#$ -q long-sl7
#$ -N dG_av_2
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l h_rt=60:00:00
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -t 1-100

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_dG_estimation.R \
	-s ${SGE_TASK_ID} \
	-m 2 \
	-d GRB2_dG_dataset_allvars

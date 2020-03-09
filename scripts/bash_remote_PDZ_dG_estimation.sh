#!/bin/bash
#$ -q long-sl7
#$ -N PDZ_dG
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l h_rt=96:00:00,h_vmem=8G,virtual_free=8G
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -t 1-10

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_PDZ_dG_estimation.R \
	-s ${SGE_TASK_ID}

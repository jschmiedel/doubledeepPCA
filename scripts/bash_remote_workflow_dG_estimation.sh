#!/bin/bash
#$ -q long-sl7
#$ -N ddPCA
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l virtual_free=55G,h_vmem=55G,h_rt=48:00:00
#$ -pe smp 15
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -t 2-5

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_workflow_dG_estimation.R -f ${SGE_TASK_ID} -l ${SGE_TASK_ID} -b 60
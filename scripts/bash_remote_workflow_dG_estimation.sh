#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N ddPCA
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l virtual_free=40G,h_vmem=40G,h_rt=24:00:00
#$ -pe smp 10
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -t 2-4

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_workflow_dG_estimation.R -f ${SGE_TASK_ID} -l ${SGE_TASK_ID} -b 10 -d GRB2

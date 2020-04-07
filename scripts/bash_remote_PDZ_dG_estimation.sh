#!/bin/bash
#$ -q long-sl7
#$ -N PDZ_dG
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -l h_rt=48:00:00,h_vmem=55G,virtual_free=55G
#$ -pe smp 15
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/doubledeepPCA/scripts/Rscript_remote_PDZ_dG_estimation.R \
		-b 0

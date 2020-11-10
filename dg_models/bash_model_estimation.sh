#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N model_est
#$ -l h_rt=6:00:00,h_vmem=4G,virtual_free=4G
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/tempura/tempura_cmdln.R \
        --dataset_folder /users/blehner/jschmiedel/doubledeepPCA/dg_models/$1 \
        --model_name $2 \
        --stage model \
        --iteration ${SGE_TASK_ID}


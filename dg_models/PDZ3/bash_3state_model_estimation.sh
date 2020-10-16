#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N PDZ3_3state_model
#$ -l h_rt=6:00:00,h_vmem=4G,virtual_free=4G
#$ -t 1:500
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

Rscript --vanilla /users/blehner/jschmiedel/tempura/tempura_cmdln.R \
        --dataset_folder /users/blehner/jschmiedel/doubledeepPCA/dg_models/PDZ3 \
        --model_name three_state \
        --stage model \
        --iteration ${SGE_TASK_ID}


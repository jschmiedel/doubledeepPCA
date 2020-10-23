#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N PDZ3_4new
#$ -l h_rt=6:00:00,h_vmem=4G,virtual_free=4G
#$ -t 1:500
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

set -eou pipefail

Rscript --vanilla /users/blehner/jschmiedel/tempura/tempura_cmdln.R \
        --dataset_folder /users/blehner/jschmiedel/doubledeepPCA/dg_models/PDZ3_newonly \
        --model_name four_state \
        --stage model \
        --iteration ${SGE_TASK_ID}


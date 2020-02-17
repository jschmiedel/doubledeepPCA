#!/bin/bash
#$ -q long-sl7,short-sl7
#$ -N Otwin_av
#$ -M jschmiedel@crg.eu
#$ -m ae
#$ -pe smp 15
#$ -l h_vmem=55G,virtual_free=55G,h_rt=24:00:00
#$ -o /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting
#$ -e /users/blehner/jschmiedel/doubledeepPCA/qsub_reporting

export JULIA_NUM_THREADS=15
# julia doubledeepPCA/scripts/JLscript_Otwinowski_prog.jl --dataset doubledeepPCA/processed_data/GRB2_dG_dataset_Otwinowski.txt --name doubledeepPCA/processed_data/GRB2_dG_Otwin

julia doubledeepPCA/scripts/JLscript_Otwinowski_prog.jl --dataset doubledeepPCA/processed_data/GRB2_dG_dataset_allvars_Otwinowski.txt --name doubledeepPCA/processed_data/GRB2_dG_allvars_Otwin

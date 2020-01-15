#!/bin/bash
#Export all environment variables
#$ -V
#Use current working directory
#$ -cwd
#Join stdout and stderr
#$ -j y
#Email address
#$ -m ae
#$ -M jschmiedel@crg.eu
#Memory
#$ -l virtual_free=50G
#Time
#$ -q long-sl7
#$ -l h_rt=24:00:00
#Parallel environment
#$ -pe smp 15

fastqFileDir="/users/project/prj004631/sequencing_data/Julia_Domingo/2019-04-18_BCDMDHANXX/"
fastqFileExtension=".fastq"
gzipped="TRUE"
experimentDesignPath="/users/blehner/jschmiedel/DMS2struct/DMS2struct_analysis/GRB2/01b-GRB2_epPCR_CYC/experimentDesign_GRB2_CYC.txt"
cutadapt5First="GGGAGGTGGAGCTAGC"
cutadapt5Second="GCGTGACATAACTAATAAGCTT"
cutadaptMinLength="90"
cutadaptErrorRate="0.2"
usearchMinQual="25"
usearchMaxee="0.5"
outputPath="/users/blehner/jschmiedel/DMS2struct/DMS2struct_analysis/GRB2/01b-GRB2_epPCR_CYC/"
projectName="DiMSum_GRB2_CYC"
startStage="1"
stopStage="0"
wildtypeSequence="ACATACGTCCAGGCCCTCTTTGACTTTGATCCCCAGGAGGATGGAGAGCTGGGCTTCCGCCGGGGAGATTTTATCCATGTCATGGATAACTCAGACCCCAACTGGTGGAAAGGAGCTTGCCACGGGCAGACCGGCATGTTTCCCCGCAATTATGTCACCCCCGTGAACTAA"
numCores="15"
fitnessMinInputCountAny="0"

#$ -N DiMSum_GRB2_CYC
#$ -t 1-1
DiMSum -i $fastqFileDir -l $fastqFileExtension -g $gzipped -e $experimentDesignPath --cutadapt5First $cutadapt5First --cutadapt5Second $cutadapt5Second  -n $cutadaptMinLength -a $cutadaptErrorRate -q $usearchMinQual -m $usearchMaxee -o $outputPath -p $projectName -s $startStage -t $stopStage -w $wildtypeSequence -c $numCores --fitnessMinInputCountAny $fitnessMinInputCountAny

#!/bin/bash
#Export all enviroment variables
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
#$ -l h_rt=40:00:00
#Parallel environment
#$ -pe smp 15

fastqFileDir="/users/project/prj004631/sequencing_data/Julia_Domingo/2019-08-12_ACDVGLANXX/"
fastqFileExtension=".fastq"
gzipped="TRUE"
experimentDesignPath="/users/blehner/jschmiedel/doubledeepPCA/dataset/GRB2/01c-GRB2_NM2_stabilityPCA/experimentDesign_GRB2_NM2_stabilityPCA.txt"
cutadapt5First="GGGAGGTGGAGCTAGC"
cutadapt5Second="GCGTGACATAACTAATAAGCTT"
cutadaptMinLength="90"
cutadaptErrorRate="0.2"
usearchMinQual="25"
usearchMaxee="0.25"
outputPath="/users/blehner/jschmiedel/doubledeepPCA/dataset/GRB2/"
projectName="01c-GRB2_NM2_stabilityPCA"
startStage="1"
stopStage="0"
wildtypeSequence="ACATACGTCCAGGCCCTCTTTGACTTTGATCCCCAGGAGGATGGAGAGCTGGGCTTCCGCCGGGGAGATTTTATCCATGTCATGGATAACTCAGACCCCAACTGGTGGAAAGGAGCTTGCCACGGGCAGACCGGCATGTTTCCCCGCAATTATGTCACCCCCGTGAACTAA"
numCores="15"
fitnessMinInputCountAny="0"
maxSubstitutions="2"

#$ -N GRB2_NM2_stabilityPCA
#$ -t 1-1
DiMSum -i $fastqFileDir -l $fastqFileExtension -g $gzipped -e $experimentDesignPath --cutadapt5First $cutadapt5First --cutadapt5Second $cutadapt5Second  -n $cutadaptMinLength -a $cutadaptErrorRate -q $usearchMinQual -m $usearchMaxee -o $outputPath -p $projectName -s $startStage -t $stopStage -w $wildtypeSequence -c $numCores --fitnessMinInputCountAny $fitnessMinInputCountAny --maxSubstitutions $maxSubstitutions



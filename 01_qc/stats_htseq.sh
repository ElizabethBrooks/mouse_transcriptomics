#!/bin/bash
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to visualize stats from htseq
# usage: sbatch stats_htseq.sh
#Submitted batch job 29907691

# required modules for servers
module load fastqc

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/counted"

# create outputs directory
outputFolder=$outputsPath"/stats_counted"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd "$outputFolder"

# name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
# add software version to output summary file
multiqc --version > $inputOutFile

# run multiqc to aggegrate the reports
multiqc -d $inputsPath"/"*"/counts.txt" -o $outputFolder -n "multiqc"

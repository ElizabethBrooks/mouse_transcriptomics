#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=mack
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to perform htseq-count counting of trimmed, aligned, then name sorted
# paired end reads
# usage: sbatch counting_htseq.sh
#Submitted batch job 

# Required modules for servers
module load htseq

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/sorted_coordinate"

#Move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="counted"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath/$outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
#Add software version to output summary file
htseq-count --version > $inputOutFile

#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsPath"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(dirname $f1)
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using htseq-count
	echo "Sample $curSampleNoPath is being counted..."
	#Use coordinate sorted flag
	#https://github.com/simon-anders/htseq/issues/37
	#--secondary-alignments ignore --supplementary-alignments ignore
	#Flag to output features in sam format
	#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
	htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i gene_id "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
	echo "Sample $curSampleNoPath has been counted!"
done

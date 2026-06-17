#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput

# script to perform htseq-count counting of trimmed, aligned, then name sorted
# paired end reads
# usage: qsub counting_htseq.sh sortedFolder
# usage Ex: qsub counting_htseq.sh sorted_coordinate
## ZQ D melanica data
## job 1945636
## EGAPx D melanica data
## job 1950551
## ZQ D melanica data
## job 2048635

#Required modules for ND CRC servers
module load bio/3.0
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2

# retrieve input folder name
sortedFolder=$1

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/genomeFeatures://g")
#genomeFile=$(grep "genomeFeatures:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/genomeFeatures://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

# setup the inputs path
inputsPath=$outputsPath"/"$sortedFolder

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
	htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
	echo "Sample $curSampleNoPath has been counted!"
done

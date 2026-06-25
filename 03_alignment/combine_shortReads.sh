#!/bin/bash
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to combine lanes of paired end reads
# Usage: sbatch combine_shortReads.sh
#Submitted batch job 

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/pairedReads://g")
# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
outputsPath=$outputsPath"/combined"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi
# move to the new directory
cd $outputsPath

# loop through all samples and combine lanes
for f1 in $readPath"/"*"_L007_R1_001.fastq.gz"; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_L007_R1_001\.fastq\.gz//')
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_L007_R1_001\.fastq\.gz//')
	# print status message
	echo "Processing $sampleTag"
	# combine lanes of forward reads
	cat $readPath"/"$sampleTag"_"*"_R1_001.fastq.gz" > $outputsPath"/"$sampleTag"_R1_001.fastq.gz"
	# combine lanes of reverse reads
	cat $readPath"/"$sampleTag"_"*"_R2_001.fastq.gz" > $outputsPath"/"$sampleTag"_R2_001.fastq.gz"
	# print status message
	echo "Processed!"
done

# print status message
echo "Analysis complete!"

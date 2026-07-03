#!/bin/bash
#SBATCH --partition=mack
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to combine lanes of paired end reads
# Usage: sbatch combine_ATAC.sh
#Submitted batch job 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# set inputs absolute path
trimmedFolder=$outputsPath"/combined_ATAC"

# make a new directory for analysis
outputsPath=$outputsPath"/combined_ATAC"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi
# move to the new directory
cd $outputsPath

# loop through all samples and combine lanes
for f1 in $trimmedFolder"/"*"_L007.R1_001.fq.gz"; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_L007.R1_001\.fq\.gz//')
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_L007.R1_001\.fq\.gz//')
	# print status message
	echo "Processing $sampleTag"
	if [[ -f $outputsPath"/"$sampleTag".R2_001.fq.gz" ]]; then
	    echo "File exists."
	else
		# combine lanes of forward reads
		cat $trimmedFolder"/"$sampleTag"_"*".R1_001.fq.gz" > $outputsPath"/"$sampleTag"_R1_001.fastq.gz"
		# combine lanes of reverse reads
		cat $trimmedFolder"/"$sampleTag"_"*".R2_001.fq.gz" > $outputsPath"/"$sampleTag"_R2_001.fastq.gz"
		# print status message
		echo "Processed!"
	fi
done

# print status message
echo "Analysis complete!"

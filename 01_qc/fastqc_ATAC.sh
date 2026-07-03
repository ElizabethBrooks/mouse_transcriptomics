#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=mack
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to perform fastqc quality control of paired end reads
# Usage: sbatch fastqc_shortReads.sh inputsType
# Usage Ex: sbatch fastqc_ATAC.sh raw
#Submitted batch job 23686947
# Usage Ex: sbatch fastqc_ATAC.sh trimmed
#Submitted batch job 23722663

# Required modules for servers
module load fastqc

# retrieve input arguments
inputsType=$1

# Retrieve paired reads absolute path for alignment
readPath=$(grep "ATACReads:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/ATACReads://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# check inputs type
if [[ $inputsType == "trimmed" ]]; then
	# set the directory for inputs
	readPath=$outputsPath"/trimmed_ATAC"
	# set the directory for analysis
	qcOut=$outputsPath"/qc_trimmed_ATAC"
elif [[ $inputsType == "raw" ]]; then
	# set the directory for analysis
	qcOut=$outputsPath"/qc_raw_ATAC"
fi

# Make a new directory for analysis
mkdir $qcOut
# Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi
# Move to the new directory
cd $qcOut

# Name version output file
versionFile=$qcOut"/software_version_summary.txt"

# Report software version
fastqc -version > $versionFile

# perform QC
fastqc $readPath"/"*".f"*"q.gz" -o $qcOut

# run multiqc
multiqc $qcOut -o $qcOut -n "multiqc"

# Print status message
echo "Analysis complete!"

#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N trimmomatic_MultiGenome_jobOutput
#$ -pe smp 8
#$ -q largemem

# Script to perform trimmomatic trimming of paired end reads
# Usage: qsub trimmomatic_shortReads.sh
## job 1843091

# Required modules for ND CRC servers
module load bio/2.0

# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
# Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/adapter://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/outputs://g")

# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir
#mkdir $outputsPath

# Make a new directory for analysis
trimOut=$outputsPath"/trimmed"
mkdir $trimOut
# Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $trimOut directory already exsists... please remove before proceeding."
	exit 1
fi
# Move to the new directory
cd $trimOut

# Name output file of inputs
versionFile=$outputsPath"/version_summary.txt"

# Add software version to outputs
echo "Trimmomatic:" >> $versionFile
trimmomatic -version >> $versionFile

# Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	# Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	# Set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# Trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	# Print status message
	echo "Processing $sampleTag"
	# phred score for trimming (Illumina 1.9)
	score=33
	# Perform adapter trimming on paired reads using 8 threads
	trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag".R1_001.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag".R2_001.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
	# Clean up
	rm -r $noPath"_R1_001_fastqc.zip"
	rm -r $noPath"_R1_001_fastqc/"
	rm -r $noPath"_R2_001_fastqc.zip"
	rm -r $noPath"_R2_001_fastqc/"
	# Print status message
	echo "Processed!"
done

# Print status message
echo "Analysis complete!"

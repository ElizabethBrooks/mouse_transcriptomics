#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=mack
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to perform samtools sorting of trimmed, then aligned
# paired end reads
# usage: sbatch sorting_ATAC_samtools.sh
#Submitted batch job 24317556
#Submitted batch job 32495027

# Required modules for servers
module load samtools

# set sorting method flags from input
methodTag="coordinate"

# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/aligned_bowtie2_ATAC"

# move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="filterd_sorted_"$methodTag"_ATAC"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
# add software version to output summary file
samtools --version > $inputOutFile

# loop through all reads and sort sam/bam files for input to samtools
for f1 in $inputsPath"/"*/; do
	# name of aligned file
	curAlignedSample="$f1"accepted_hits.bam
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	# create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	# status message
	echo "Sample $curSampleNoPath is being processed..."
	# remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates
	samtools view -@ 8 -F 524 -f 2 -u "$curAlignedSample" | samtools sort -@ 8 -n -o "$outputFolder"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam
	# update paired-end flags for singletons, and remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
	samtools fixmate -@ 8 -m -r "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
	# remove mitochondiral sequences and obtain position sorted BAM
	samtools view -@ 8 -F 1804 -f 2 -u "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam | grep -v -E '\tMT\t' | samtools sort -@ 8 -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam
	# clean up
	rm "$outputFolder"/"$curSampleNoPath"/sortedName.bam
	rm "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
	# remove duplicate reads
	samtools markdup -@ 8 -r "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam "$outputFolder"/"$curSampleNoPath"/markedDups.bam
	# index bam files
	samtools index -@ 8 "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
	samtools index -@ 8 "$outputFolder"/"$curSampleNoPath"/markedDups.bam
	# status message
	echo "Sample $curSampleNoPath has been processed!"
done

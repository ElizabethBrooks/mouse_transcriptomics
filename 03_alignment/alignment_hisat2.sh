#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=mack
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to perform hisat2 alignment of trimmed
# paired end reads
# Note that a hisat2 genome refernce build folder needs to be generated first
# usage: sbatch alignment_hisat2.sh
#Submitted batch job 23685130

# Required modules for servers
module load hisat2

#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../inputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/pairedReads://g")
# retrieve input trimmed reads path
inputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# set inputs absolute path
trimmedFolder=$inputsPath"/combined"

# move to outputs directory
cd "$outputsPath"

# set output directory name
outputFolder=$outputsPath"/aligned"
# create output directory
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
#Add software version to output summary file
hisat2 --version > $inputOutFile
samtools --version >> $inputOutFile

#Build output directory for Hisat reference
buildOut="$outputsPath"/"reference_hisat2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//' | sed 's/\.fna//' | sed 's/\.fa//')

#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in $trimmedFolder"/"*_R1_001.fastq.gz; do
	# status message
	echo "Processing file $f1 ..."
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R1_001\.fastq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/_R1_001\.fastq\.gz//')
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Run hisat2 with default settings
	echo "Sample $curSampleNoPath is being aligned and converted..."
	hisat2 -p 8 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample".R2_001.fastq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam \
	--un-conc-gz "$outputFolder"/"$curSampleNoPath"/un_conc.fq.gz --al-conc-gz "$outputFolder"/"$curSampleNoPath"/al_conc.fq.gz --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 8 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
	#Remove the now converted .sam file
	rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	# status message
	echo "Sample $curSampleNoPath has been aligned and converted!"
done

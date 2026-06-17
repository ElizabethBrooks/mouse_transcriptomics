#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_test_jobOutput
#$ -pe smp 4

# Script to perform hisat2 alignment of trimmed
# paired end reads
# Note that a hisat2 genome refernce build folder needs to be generated first
# usage: qsub alignment_hisat2_test.sh
# run 1
# test
## job 1928484
# EGAPx test
## job 1980157
# run 2
# test 
## job 2050275
# EGAPx test
## job 2050277
# un-conc and al-conc
# test
## job 2071558
# EGAPx test
## job 2071554

#Required modules for ND CRC servers
module load bio/2.0
#module load bio/hisat2/2.1.0

#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
#buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_EGAPx_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_D_pulex.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_EGAPx_D_pulex.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")
# retrieve input trimmed reads path
inputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/outputs://g")

# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir
inputsPath=$inputsPath"/"$projectDir

# set inputs absolute path
trimmedFolder=$inputsPath"/trimmed"

# move to outputs directory
cd "$outputsPath"

# set output directory name
outputFolder=$outputsPath"/aligned_conc"
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

#Build output directory for Hisat reference
buildOut="$outputsPath"/"reference_hisat2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//' | sed 's/\.fna//' | sed 's/\.fa//')

#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in $trimmedFolder"/"*.R1_001.fq.gz; do
	# status message
	echo "Processing file $f1 ..."
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/\.R1_001\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/\.R1_001\.fq\.gz//')
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Run hisat2 with default settings
	echo "Sample $curSampleNoPath is being aligned and converted..."
	hisat2 -p 4 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample".R2_001.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam \
	--un-conc-gz "$outputFolder"/"$curSampleNoPath"/un_conc.fq.gz --al-conc-gz "$outputFolder"/"$curSampleNoPath"/al_conc.fq.gz --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
	#Remove the now converted .sam file
	rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	# status message
	echo "Sample $curSampleNoPath has been aligned and converted!"
done

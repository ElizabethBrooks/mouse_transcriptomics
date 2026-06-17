#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sorting_samtools_jobOutput
#$ -pe smp 4

# Script to perform samtools sorting of trimmed, then aligned
# paired end reads
# usage: qsub sorting_samtools.sh sortingMethod alignedFolder optionalAssembledFolder
# usage Ex: qsub sorting_samtools.sh name
## ZQ D melanica data
## job 1931639
# usage Ex: qsub sorting_samtools.sh coordinate
## ZQ D melanica data
## job 1943087
# usage Ex: qsub sorting_samtools.sh coordinate
## EGAPx D melanica data
## job 1949843

#Required modules for ND CRC servers
module load bio

#Retrieve sorting method flags from input
if [[ "$1" == "name" || "$1" == "Name" || "$1" == "n" || "$1" == "N" ]]; then
	#Name sorted flag with num threads flag
	flags="-@ 4 -n"
	methodTag="name"
elif [[ "$1" == "coordinate" || "$1" == "Coordinate" || "$1" == "c" || "$1" == "C" ]]; then
	#Coordinate sorted with num threads flag
	flags="-@ 4"
	methodTag="coordinate"
else
	#Report error with input flag
	echo "ERROR: a flag for sorting method (name or coordinate) is expected... exiting"
	exit 1
fi

# Retrieve analysis outputs absolute path
#outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

# setup the inputs path
inputsPath=$outputsPath"/aligned"

#Move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="sorted_"$methodTag
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Add software version to output summary file
samtools --version > $inputOutFile

#Loop through all reads and sort sam/bam files for input to samtools
for f1 in $inputsPath"/"*/; do
	#Name of aligned file
	curAlignedSample="$f1"accepted_hits.bam
	#Trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Run samtools to prepare mapped reads for sorting by name
	#using 8 threads
	echo "Sample $curSampleNoPath is being name sorted..."
	samtools sort -@ 4 -n -o "$outputFolder"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample"
	echo "Sample $curSampleNoPath has been name sorted!"
	#Determine which sorting method is to be performed
	if [[ "$methodTag" == "coordinate" ]]; then
		#Run fixmate -m to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate -m "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		#Clean up
		rm "$outputFolder"/"$curSampleNoPath"/sortedName.bam
		#Run samtools to prepare mapped reads for sorting by coordinate
		#using 8 threads
		echo "Sample $curSampleNoPath is being sorted..."
		samtools sort "$flags" -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath has been sorted!"
		rm "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		#Remove duplicate reads
		samtools markdup -r "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam "$outputFolder"/"$curSampleNoPath"/markedDups.bam
		# index bam files
		samtools index -@ 4 "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
		samtools index -@ 4 "$outputFolder"/"$curSampleNoPath"/markedDups.bam
	else
		#Run fixmate -m to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate -m "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		rm "$outputFolder"/"$curSampleNoPath"/sortedName.bam
		#Remove duplicate reads
		#samtools markdup -r "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam "$outputFolder"/"$curSampleNoPath"/markedDups.bam
	fi
done
#Copy previous summaries
#cp "$inputsPath"/"$3"/*.txt "$outputFolder"

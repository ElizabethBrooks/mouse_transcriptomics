#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N building_hisat2_test_jobOutput
#$ -pe smp 4

#Script to generate a hisat2 genome refernce build folder
# usage: qsub building_hisat2_test.sh
# test
## job 1928389
# EGAPx
## job 1980152

#Required modules for ND CRC servers
module load bio/2.0
#module load bio/hisat2/2.1.0

#Retrieve genome reference absolute path for alignment
#buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_EGAPx_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve build outputs absolute path
#outputsPath=$(grep "outputs:" ../inputData/shortReads/inputPaths_D_pulex.txt | tr -d " " | sed "s/outputs://g")
outputsPath=$(grep "outputs:" ../inputData/shortReads/inputPaths_EGAPx_D_pulex.txt | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")

# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Move to outputs directory
cd "$outputsPath"

#Create output directory
outputFolder=$outputsPath"/reference_hisat2_build"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
#Add software version to output summary file
hisat2-build --version > $inputOutFile

# check if the build already exists
if [ $? -eq 0 ]; then
	#Trim file path from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNewPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g' | sed 's/\.fna/\.fa/g')
	#Copy genome build fasta file to hisat2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNewPath"
	#Trim file extension
	buildFileNoPath=$(echo $buildFileNewPath | sed 's/\.fa//g')
	#Begin hisat2 build
	echo "Beginning hisat2 build... "
	hisat2-build -p 4 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath"
	echo "hisat2 build complete!"
else
	echo "Build folder reference_hisat2_build already exists, skipping building..."
fi

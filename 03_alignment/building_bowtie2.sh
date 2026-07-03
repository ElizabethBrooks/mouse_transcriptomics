#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to generate a hisat2 genome refernce build folder
# usage: sbatch building_bowtie2.sh
#Submitted batch job 23725016

# Required modules for servers
module load bowtie2

#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../inputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

#Move to outputs directory
cd "$outputsPath"

#Create output directory
outputFolder=$outputsPath"/reference_bowtie2_build"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
#Add software version to output summary file
bowtie2-build --version > $inputOutFile

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
	echo "Beginning bowtie2 build... "
	bowtie2-build --threads 8 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath"
	echo "bowtie2 build complete!"
else
	echo "Build folder reference_hisat2_build already exists, skipping building..."
fi

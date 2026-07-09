#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=mack
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to generate alignment stats using samtools
# paired end reads
# usage: sbatch stats_samtools.sh inputsType
# usage ex: sbatch stats_ATAC_samtools.sh hisat2
#Submitted batch job 
# usage ex: sbatch stats_ATAC_samtools.sh bowtie2
#Submitted batch job 

# required modules for servers
module load samtools

# retrieve input arguments
inputsType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# check inputs type
if [[ $inputsType == "hisat2" ]]; then
	# setup the inputs path
	inputsPath=$outputsPath"/aligned_ATAC"
	# set the directory for analysis
	outputFolder=$outputsPath"/stats_aligned_ATAC"
elif [[ $inputsType == "bowtie2" ]]; then
	# setup the inputs path
	inputsPath=$outputsPath"/aligned_bowtie2_ATAC"
	# set the directory for analysis
	outputFolder=$outputsPath"/stats_aligned_bowtie2_ATAC"
fi

# create outputs directory
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd "$outputFolder"

# name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
# add software version to output summary file
samtools --version > $inputOutFile

# loop through all reads and sort bam files for input to samtools
for f1 in $inputsPath"/"*/; do
	# name of aligned file
	curAlignedSample="$f1"accepted_hits.bam
	#Trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	# run samtolls stats
	samtools stats -@ 8 $curAlignedSample > $outputFolder"/"$curSampleNoPath".stats"
done

# run multiqc to aggegrate the reports
multiqc $outputFolder -o $outputFolder -n "multiqc"

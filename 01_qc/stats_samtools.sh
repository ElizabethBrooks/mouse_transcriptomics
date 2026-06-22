#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to generate alignment stats using samtools
# paired end reads
# usage: sbatch stats_samtools.sh
# usage ex: sbatch stats_samtools.sh
#Submitted batch job 

# required modules for servers
module load samtools

# retrieve sorting method flags from input
if [[ "$1" == "name" || "$1" == "Name" || "$1" == "n" || "$1" == "N" ]]; then
	# name sorted flag with num threads flag
	flags="-@ 4 -n"
	methodTag="name"
elif [[ "$1" == "coordinate" || "$1" == "Coordinate" || "$1" == "c" || "$1" == "C" ]]; then
	# coordinate sorted with num threads flag
	flags="-@ 4"
	methodTag="coordinate"
else
	# report error with input flag
	echo "ERROR: a flag for sorting method (name or coordinate) is expected... exiting"
	exit 1
fi

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/aligned"

# move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="aligned_stats"
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

# loop through all reads and sort bam files for input to samtools
for f1 in $inputsPath"/"*/; do
	# name of aligned file
	curAlignedSample="$f1"accepted_hits.bam
	# run samtolls stats
	samtools stats -@ 8 $curAlignedSample
done

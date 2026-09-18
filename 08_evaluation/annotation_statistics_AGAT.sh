#!/bin/bash
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to keep only the longest isoforms in the input gff
# usage: sbatch annotation_statistics_AGAT.sh

# load software
conda activate my_agat

# retrieve genome reference and features absolute paths
genomeRef=$(grep "genomeReference:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="AGAT"
mkdir "$outputFolder"

# create outputs directory
mkdir $outputsPath"/AGAT"

# status message
echo "Beginning analysis of $speciesName..."

# pre clean
rm $outputsPath"/AGAT/annotation_stats.txt"

# extract annotation statistics
agat_sp_statistics.pl --gff $genomeFeat -g $genomeRef -o $outputsPath"/AGAT/annotation_stats.txt"

# status message
echo "Analysis complete!"

#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to keep only the longest transcripts
# usage: sbatch extract_longest_transcripts_AGAT.sh

# load software
conda activate my_agat

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")

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
rm $outputsPath"/AGAT/longest_mRNA.fa"

# extract longest transcripts
agat_sp_extract_sequences.pl -gff $outputsPath"/AGAT/output_longest.gff" -f $genomeFile -t mRNA -o $outputsPath"/AGAT/longest_mRNA.fa" --thread 8

# status message
echo "Analysis complete!"

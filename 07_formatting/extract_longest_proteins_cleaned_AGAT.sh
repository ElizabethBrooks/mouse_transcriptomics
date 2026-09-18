#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to keep only the longest proteins
# usage: sbatch extract_longest_proteins_cleaned_AGAT.sh

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
rm $outputsPath"/AGAT/longest_protein_cleaned.fa"

# extract longest proteins
agat_sp_extract_sequences.pl -g $outputsPath"/AGAT/output_longest.gff" -f $genomeFile -p -o $outputsPath"/AGAT/longest_protein_cleaned.fa" --clean_final_stop --clean_internal_stop --thread 8

# status message
echo "Analysis complete!"

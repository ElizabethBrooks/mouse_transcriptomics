#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to keep only the longest isoforms in the input gff
# usage: sbatch keep_longest_isoforms_AGAT.sh
## Submitted batch job 29903589

# load software
conda activate my_agat

# retrieve genome features absolute path
genomeFile=$(grep "genomeFeatures:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="AGAT"
mkdir "$outputFolder"

# create outputs directory
mkdir $outputsPath"/AGAT"

# name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
# add software version to output summary file
agat --version > $inputOutFile

# status message
echo "Beginning analysis..."

# pre clean
rm $outputsPath"/AGAT/output_longest.gff"

# extract longest isoforms
agat_sp_keep_longest_isoform.pl -f $genomeFile -o $outputsPath"/AGAT/output_longest.gff" --thread 8

# status message
echo "Analysis complete!"

#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# script to extract proteins
# usage: sbatch extract_proteins_cleaned_gffread.sh

# retrieve software path
softwarePath="/kuhpc/scratch/mack/e959b751/software/gffread"

# retrieve genome reference features absolute path
genomeRef=$(grep "genomeReference:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")

# retrieve analysis outputs absolute path
outputsPath=$(dirname $genomeFile)

# move to software directory
cd $softwarePath

# status message
echo "Beginning analysis..."

# extract proteins
./gffread -y $outputsPath"/Mus_musculus.GRCm39.proteins.fa" -g $genomeRef $genomeFeat -S -V

# status message
echo "Analysis complete!"

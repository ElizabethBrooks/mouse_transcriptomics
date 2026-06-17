#!/bin/bash

# script to generate a merge gene counts file
# usage: bash generateMergedTable.sh

# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/counted"

# initialize the merged counts file
echo "gene" > $inputsPath"/counts_merged.tmp.csv"
firstFile=$(ls -1 $inputsPath | head -n 1)
cat $firstFile"/counts.txt" | cut -f1 >> $inputsPath"/counts_merged.tmp.csv"

# merge counts for each sample
for i in $inputsPath"/"*"/"; do 
	# clean up sample name
	newName=$(basename $i | sed "s/_S.*_L004//g")
	# status message
	echo "Processing sample $newName..."
	# add sample name to the sample outputs
	echo $newName > $inputsPath"/counts_sample.tmp.csv"
	# add counts to the sample outputs
	cat $i/counts.txt | cut -f2 >> $inputsPath"/counts_sample.tmp.csv"
	# merge sample counts
	paste -d, $inputsPath"/counts_merged.tmp.csv" $inputsPath"/counts_sample.tmp.csv" > $inputsPath"/counts_merged.csv"
	# update the merged counts data
	cat $inputsPath"/counts_merged.csv" > $inputsPath"/counts_merged.tmp.csv"
done

# clean up
rm $inputsPath"/counts_sample.tmp.csv"
rm $inputsPath"/counts_merged.tmp.csv"

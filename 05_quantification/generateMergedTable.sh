#!/bin/bash

# script to generate a merge gene counts file
# usage: bash generateMergedTable.sh

# set the inputs path
#inputsPath="/scratch365/ebrooks5/D_melanica_UV_exposure/short_read_data_processed_test/Pfrender_MP-3533_250512_CMG/counted"
#inputsPath="/scratch365/ebrooks5/D_melanica_UV_exposure/short_read_data_processed_EGAPx/Pfrender_MP-3533_250512_CMG/counted"
inputsPath="/scratch365/ebrooks5/D_melanica_UV_exposure/short_read_data_processed_EGAPx_test/Pfrender_MP-3533_250512_CMG/counted"

# initialize the merged counts file
echo "gene" > $inputsPath"/counts_merged.tmp.csv"
cat $inputsPath"/MUV10_S10_L004/counts.txt" | cut -f1 >> $inputsPath"/counts_merged.tmp.csv"

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

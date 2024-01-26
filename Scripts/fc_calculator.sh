#!/usr/bin/bash

IFS=$'\t'

#Remaking a list of all the possible conditions for a sample group, including any future additions
sample_types=$(cut -f2 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "SampleType")
time_types=$(cut -f4 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "Time")
treatment_types=$(cut -f5 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "Treatment")


#User-generated group-wise comparisons to display fold changes across gene expression counts
mkdir fold_changes

echo -e "What group-wise comparison would you like to make for gene expression fold change? Please input 6 condition names of the two groups in the same format as displayed below, each separated by a comma, and in this order: two for sample type, two for treatment type, and two for time type. Ensure the first group to compare is described in the first of each condition inputted, while the second group to compare is described in the second of each condition inputted. Note that the fold change will be calculated by dividing the first set by the second set.  For example, if you are comparing fold change across 24h vs 48h of Induced of WT samples, you would input 'WT,WT,Induced,Induced,48,24'"
echo -e "These are the available sample types: \n${sample_types}"
echo -e "These are the available treatment types: \n${treatment_types}"
echo -e "These are the available time types: \n${time_types}"

read -p "Please type your 6 condition names:" input

IFS=$','

#Assigns the inputs into separate variables
read -r sample_type_1 sample_type_2 treatment_type_1 treatment_type_2 time_type_1 time_type_2 <<< "${input}"

IFS=$'\t'

#Checking the sample groups files have passed quality check and average counts of gene expression have been generated
if [[ -f ./counts/counts_${sample_type_1}_${treatment_type_1}_${time_type_1}_averaged && -f ./counts/counts_${sample_type_2}_${treatment_type_2}_${time_type_2}_averaged ]]
	then
		#Calculates the fold change between the two sample groups inputted by the user and stores this in a new file, alongside the gene name and gene description. If a gene expression count is 0, a small value is added to both counts on the numerator and denominator to ensure that a fold change can be calculated.
		paste -d '\t' ./counts/counts_${sample_type_1}_${treatment_type_1}_${time_type_1}_averaged ./counts/counts_${sample_type_2}_${treatment_type_2}_${time_type_2}_averaged | awk 'BEGIN{FS="\t"; OFS="\t";}{if ($1 == 0 || $4 == 0){fc=($1 + 0.0001)/($4 + 0.0001) ;} else {fc= $1/$4 ;} print fc, $2, $3;}' | sort -k1,1 -r > ./fold_changes/${sample_type_1}_${treatment_type_1}_${time_type_1}_vs_${sample_type_2}_${treatment_type_2}_${time_type_2}_fold_change
		
		echo "Below is the first few lines of the generated output files. Your file with the calculated fold changes has been saved in ./fold_changes/${sample_type_1}_${treatment_type_1}_${time_type_1}_vs_${sample_type_2}_${treatment_type_2}_${time_type_2}_fold_change"
		echo "$(head ./fold_changes/${sample_type_1}_${treatment_type_1}_${time_type_1}_vs_${sample_type_2}_${treatment_type_2}_${time_type_2}_fold_change)"

	else
		echo "Please double check your spelling or format for the conditions. It seems the gene expression counts file for your query is not available, it has likely failed quality check."

fi

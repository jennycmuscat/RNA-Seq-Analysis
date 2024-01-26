#!/bin/bash

#Copying the fastq files to the user's current local directory
scp -r /localdisk/data/BPSM/ICA1/fastq .

#Change the permission so that the following commands can access the fastq files
chmod 700 fastq 

#Creating a directory to save the output of quality checks
mkdir ./fastqc_output
touch ./fastqc_output/accepted_samples
chmod 700 ./fastqc_output/accepted_samples

unset IFS
IFS=$'\t'


#Running quality checks for each of the fastq files logged on Tco2.fqfiles. This is done for both end1 and end2
#Asks the user to input the maximum threshold of fails from fastqc
read -p "What is the maximum threshold of fails out of 10 for the fastq file to pass quality check?" threshold
while read sample_name sample_type replicate time treatment end1 end2 
do
	if [[ ${sample_name} != "SampleName" ]]
	then
		for num in {1..2}
		do
			#Creates the file name of the sample from the sample name. For example Tco999 becomes Tco-999
			name="$(echo "${sample_name}" | cut -c 1-3)-$(echo "${sample_name}" | cut -c 4-)"
			fastqc -o ./fastqc_output ./fastq/${name}_${num}.fq.gz
			unzip -o -q -d ./fastqc_output ./fastqc_output/${name}_${num}_fastqc.zip
			
			
			#Checks the summary.txt output of fastqc to determine the total number of quality check fails each sample got
			fails_count=0
			while read status check file
			do
				if [[ ${status} == "FAIL" ]]
				then
					((fails_count++))
				fi
			done < ./fastqc_output/${name}_${num}_fastqc/summary.txt


			if [[ ${fails_count} > ${threshold} ]]
			then
				#If even the first pair end read fails quality check, neither of the pair ends will process for analysis
				echo "Sample ${sample_name} has failed quality check and will not be analysed."
				break
				#Ensuring that if one of the pair ends fails qc, then neither files for end #1 or end #2 will proceed to analysis
			elif [[ ${num} == 2 ]]
			then
				#Adds the name of the samples that pass quality check to a file for future use
				echo -e "${name}" >> ./fastqc_output/accepted_samples
			fi
		done 
	else
		continue
	fi
done < ./fastq/Tco2.fqfiles



#Copy the reference genome to local directory
scp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .

#Index the reference genome into .bt2 format saved in the Tcongo_genome directory
bowtie2-build ./Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz ./Tcongo_genome/Tcongo_genome_bt2

#Making a directory for the alignment read pairs
mkdir alignment_reads
chmod 700 ./fastqc_output/accepted_samples


#Loops through the samples that passed quality check and aligns both end pairs to the reference, outputting a bam file in the directory "algnment_reads" alongside their index files (".bam.bai")
while read file
do
	bowtie2 -x ./Tcongo_genome/Tcongo_genome_bt2 -1 ./fastq/${file}_1.fq.gz -2 ./fastq/${file}_2.fq.gz | samtools view -bS | samtools sort -o ./alignment_reads/${file}_sorted.bam
	samtools index ./alignment_reads/${file}_sorted.bam
done < ./fastqc_output/accepted_samples



#Making a directory to save the counts data and reference index in
mkdir counts
scp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ./counts

IFS=$'\t'

#Making a list of all the possible conditions for a sample group, including any future additions
sample_types=$(cut -f2 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "SampleType")
time_types=$(cut -f4 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "Time")
treatment_types=$(cut -f5 ./fastq/Tco2.fqfiles | sort | uniq | grep -v "Treatment")

IFS=$'\n'


#Generating counts data and saving these in separate txt files, each containing one sample group with all replicates
echo "${sample_types}" | while read -r types
do
	echo "${time_types}" | while read -r time
	do
		echo "${treatment_types}" | while read -r treatment
		do
			IFS=$'\t'
			files=""
			
	
			#Loops through all the samples logged on Tco2.fqfiles and identifies the replicates of each sample group type
			while read sample_name sample_type replicate times treatments end1 end2
			do
				if [[ ${sample_type} == ${types} && ${times} == ${time} && ${treatments} == ${treatment} ]]
				then
					#Creates the file name from the name given on Tco2.fqfiles. For example Tco999 becomes Tco-999
					path="$(echo "${sample_name}" | cut -c 1-3)-$(echo "${sample_name}" | cut -c 4-)"
					#Checks that this sample name is part of the list of samples that passed quality check
					if grep -Fxq "${path}" ./fastqc_output/accepted_samples
					then
						#Concatenating all the replicates file paths of alignment read pairs of one sample group into one variable, separated by a space
						IFS=$','
						files+="./alignment_reads/${path}_sorted.bam,"
			
						#Splits the variable containing all file paths of the alignement read pairs into separate variables (up to a limit of 12 replicates)
						read -r a b c d e f g h i j k l <<< ${files}
						
						#Generates the actual counts files, with a file per sample group containing the counts of all replicates in a separate column
						bedtools multicov -bams ${a} ${b} ${c} ${d} ${e} ${f} ${g} ${h} ${i} ${j} ${k} ${l} -bed ./counts/TriTrypDB-46_TcongolenseIL3000_2019.bed > ./counts/counts_${types}_${treatment}_${time}
						IFS=$'\t'

					fi
				fi
			done < ./fastq/Tco2.fqfiles


			IFS=$'\t'

			#Calculates the number of repliactes present in each sample group file by coutning the columns other than the standard 5 gene columns as well as disregards columns with null content in them
			if [[ -f ./counts/counts_${types}_${treatment}_${time} ]]
                        then
				row=$(head -n1 ./counts/counts_${types}_${treatment}_${time})
				replicates=$(awk 'BEGIN{FS="\t"; OFS="\t"}{count=0; for (i=1; i<=NF; i++) if (length($i) != ""){count++;} print count - 5;}'  <<< "${row}")
		
				#Creates a file for each sample group containing the gene name, gene description and average counts across replicates
				awk 'BEGIN{FS="\t"; OFS="\t";}{sum=0; i=1; while (i <= re){sum+=$(5+i); i++;} print sum/re, $4, $5;}' re=${replicates} ./counts/counts_${types}_${treatment}_${time} > ./counts/counts_${types}_${treatment}_${time}_averaged
			fi
		done
	done
done


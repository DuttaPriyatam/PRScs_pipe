#!/bin/bash

helpFunction()
{
	echo ""
	echo "Usage: $0 -i input -o output"
	echo -e "\t-i Input file with the path to directory "
	echo -e "\t-o Output file with phenotype Abbreviation(e.g., /path/to/dir/SCZ for Schizophrenia)"
	exit 1 #Exit script after printing help
}

while getopts "i:o:" opt
do
	case "$opt" in
		i ) input="$OPTARG" ;;
		o ) output="$OPTARG" ;;
		? ) helpFunction ;; #Print helpFunction in case parameter is non-existent
	esac
done

# Print helpFunction in case parameters are emplty
if [ -z "$input" ] || [ -z "$output" ]
then
	echo "Some or all of the parameters are empty"
	helpFunction
fi

# Begin Scripts in case all parameters are correct
echo "Input: $input"
echo "Output: $output"

# Removing meta_information from the base sum stat file
grep -v "#" $input > "$output"_base.vcf.tsv
echo "Total SNPs: $(sed '1d' ${output}_base.vcf.tsv | wc -l)"

# Here A1 - effect allele(reference allele) and A2 is the altenate or non-effect allele.

# The file has already been filtered for high imputation quality(INFO >0.8), and MAF > 0.01.
# MAF calculated using the following formula: (FCAS*NCAS)+(FCON*$NCON))/($NCON+$NCAS)
#awk -F'\t' 'NR==1 || ($8 > 0.8 && $6 > 0.01 && $7 > 0.01) {print}' "$output"_base.vcf.tsv > "$output"_filtered.vcf.tsv
awk -F'\t' 'NR==1 || ($8 > 0.8 && (($6*$12)+($7*$13))/($12+$13) > 0.01) {print}' "$output"_base.vcf.tsv > "$output"_filtered.vcf.tsv
echo "SNPs after filtering for Imputation and MAF: $( sed '1d' ${output}_filtered.vcf.tsv | wc -l)"

# Remove duplicates - The current file doesn't have any duplicate
awk -F'\t' '{seen[$2]++; if(seen[$2]==1){ print}}' "$output"_filtered.vcf.tsv > "$output"_nodup.vcf.tsv
echo "SNPs after removing duplicated SNPs: $(sed '1d' ${output}_nodup.vcf.tsv | wc -l)"

# Removing ambiguous SNPS
awk -F'\t' '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' "$output"_nodup.vcf.tsv  > "$output"_QC.vcf.tsv
echo "SNPs after correcting for ambiguous SNPs/Final variants after QC: $(sed '1d' ${output}_QC.vcf.tsv | wc -l)"

# Updating the file for PRScs
sed -i -e "s/CHROM/CHR/g" -e "s/POS/BP/g" -e "s/ID/SNP/g" "$output"_QC.vcf.tsv | cut -f2,4-5,9-10 > "$output"_for_PRScs.vcf.tsv
echo "${output}_for_PRScs.vcf.tsv created for PRScs calculation"

# Removing intermediate files
rm "$output"_filtered.vcf.tsv
rm "$output"_nodup.vcf.tsv
rm "$output"_base.vcf.tsv

echo "Intermediate files removed"

# Extract rsids from the QCed sum stats and create the rsid list file
sed '1d' "$output"_QC.vcf.tsv | cut -f2 > "$output"_rsid.txt
echo "Base rsid extracted to ${output}_rsid.txt"

# Creating a csv file for target QC
sed -e "s/CHROM/CHR/g" -e "s/\t/,/g" "$output"_QC.vcf.tsv > "$output"_QC.csv

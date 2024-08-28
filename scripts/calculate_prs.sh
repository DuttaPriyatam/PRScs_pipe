#!/bin/bash

helpFunction()
{
	echo ""
	echo "Usage: $0 -b base -t target -r ref -p phe -n num -o output"
	echo -e "\t-b Path to the base summ stat file"
	echo -e "\t-t Path to target bim file"
	echo -e "\t-r Path to the LD Reference panel"
	echo -e "\t-p Phenotype abbreviation(e.g., SCZ for Schizophrenia)"
	echo -e "\t-n Total GWAS participants"
	echo -e "\t-o Output directory folder"
	exit 1 #Exit script after printing help
}

while getopts "b:t:r:p:n:o:" opt
do
	case "$opt" in
		b ) base="$OPTARG" ;;
		t ) target="$OPTARG" ;;
		r ) ref="$OPTARG" ;;
		p ) phe="$OPTARG" ;;
		n ) num="$OPTARG" ;;
		o ) output="$OPTARG" ;;
		? ) helpFunction ;; #Print helpFunction in case parameter is non-existent
	esac
done

# Print helpFunction in case parameters are emplty
if [ -z "$base" ] || [ -z "$target" ] || [ -z "$ref" ] || [ -z "$phe" ] || [ -z "$num" ] || [ -z "$output" ]
then
	echo "Some or all of the parameters are empty"
	helpFunction
fi

declare -a phi_vals=("1" "1e-2" "1e-4" "1e-6")

# PRScs auto
python /home/pdutta/scratch/apps/PRS-CS/PRScs/PRScs.py --ref_dir="$ref" --bim_prefix="$target" --sst_file="$base" --n_gwas="$num" --out_dir="$output"auto

# Concatenating the output files
cat "$output"auto_pst* > "$output"auto_prscs.txt
# Removing intermediate files
rm "$output"/auto_pst*
# Computing PRScs
plink2 --bfile "$target" --memory 16000 --score "$output"auto_prscs.txt 2 4 6 --out "$output"auto

# # PRScs for different phi values
# for val in "${phi_vals[@]}"
# do
# 	python /home/pdutta/scratch/apps/PRS-CS/PRScs/PRScs.py --ref_dir="$ref" --bim_prefix="$target" --sst_file="$base" --n_gwas="$num" --phi="$val" --out_dir="$output"
	
# 	# Concatenating the output files
# 	cat "$output"*phi"$val"* > "$output"phi_"$val"_prscs.txt
	
# 	# Removing intermediate files
# 	rm "$output"/_pst*

# 	# Computing PRScs
# 	plink2 --bfile "$target" --memory 14000 --score "$output"phi_"$val"_prscs.txt 2 4 6 --out "$output"_phi_"$val"

# done

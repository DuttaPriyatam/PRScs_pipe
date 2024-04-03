#!/bin/bash

helpFunction()
{
	echo ""
	echo "Usage: $0 -imp input -base base -phe phenotype -out output"
	echo -e "\t-imp Path to imputation files"
	echo -e "\t-base Path to QC'd base dataset files"
	echo -e "\t-phe Phenotype abbreviation(e.g., SCZ for Schizophrenia)"
	echo -e "\t-out Output directory folder"
	exit 1 #Exit script after printing help
}

while getopts "imp:base:phe:out:" opt
do
	case "$opt" in
		imp ) imp_file_dir="$OPTARG" ;;
		base ) base="$OPTARG" ;;
		phe ) phenotype="$OPTARG" ;;
		out ) output="$OPTARG" ;;
		? ) helpFunction ;; #Print helpFunction in case parameter is non-existent
	esac
done

# Print helpFunction in case parameters are emplty
if [ -z "$imp_file_dir" ] || [ -z "$base" ] || [ -z "$phenotype" ] || [ -z "$output" ]
then
	echo "Some or all of the parameters are empty"
	helpFunction
fi

#Imputation file path
imp_file_dir="/home/pdutta/scratch/data/imputation"

#Output folder
output="/home/pdutta/scratch/data/SCZ/target/filter_1"

#Base folder
base="/home/pdutta/scratch/data/SCZ/base"

# Script folder
scripts="/home/pdutta/scratch/scripts"

for i in {20..22}
do
	bgenix -g $imp_file_dir/ukb22828_c${i}_b0_v3.bgen -incl-rsids $base/rsid.txt > $output/chr${i}.bgen 
done

#Concatenating the extreacted bgen files
cat-bgen -g $output/*.bgen -og $output/initial_chr.bgen -clobber

# Writing the index
bgenix -g $output/initial_chr.bgen -index -clobber

for i in {20..22}
do
	rm $output/chr${i}.bgen
done

sed -e "s/\t/,/g" $base/SCZ_QC.vcf.tsv > $base/SCZ_QC.csv

# Import the betas into the sqlite database as a table called Betas
sqlite3 $output/initial_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," $output/initial_chr.bgen.bgi ".import $base/SCZ_QC.csv Betas"

sqlite3 $output/initial_chr.bgen.bgi "DROP TABLE IF EXISTS Joined"
# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header -csv $output/initial_chr.bgen.bgi \
"CREATE TABLE Joined AS
 SELECT Variant.*, Betas.CHR, Betas.BETA FROM Variant INNER JOIN Betas
   ON Variant.chromosome = printf('%02d', Betas.CHR)
   AND Variant.position = Betas.POS
   AND Variant.allele1 = Betas.A2
   AND Variant.allele2 = Betas.A1
 UNION
 SELECT Variant.*, Betas.CHROM, -Betas.BETA FROM Variant INNER JOIN Betas
   ON Variant.chromosome = printf('%02d', Betas.CHR)
   AND Variant.position = Betas.POS
   AND Variant.allele1 = Betas.A1 AND
   Variant.allele2 = Betas.A2;"

# Filter the .bgen file to include only the alleles specified in the Betas for each SNP
bgenix -g $output/initial_chr.bgen -table Joined  > $output/single_allelic.bgen

# And produce an index file for the new .bgen
bgenix -g $output/single_allelic.bgen -index


#--set-all-var-ids @:#_\$r\$a \
# Converting file to bed, bim and fam format
plink2 \
	--bgen $output/single_allelic.bgen ref-first \
 	--sample $imp_file_dir/ukb22828_c22_b0_v3.sample \
 	--memory 16000 \
 	--threads 8 \
	--freq \
 	--make-bed \
 	--out $output/raw

#Marking ambigous SNPs from target data (I don't think step is relevant)
awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' $output/raw.afreq > $output/exclrsIDs_ambiguous.txt


# Standard GWAS QC
plink2 --bfile $output/raw \
  --memory 16000 \
 	--threads 8 \
	--exclude $output/exclrsIDs_ambiguous.txt \
	--maf 0.01 \
	--hwe 1e-6 \
	--geno 0.01 \
	--mind 0.01 \
  --rm-dup "exclude-mismatch" \
	--write-snplist \
	--make-just-fam \
	--out $output/raw.QC

# # Pruning to remove highly correlated SNPs
plink2 \
    --bfile $output/raw \
    --memory 16000 \
 	  --threads 8 \
    --keep $output/raw.QC.fam \
    --extract $output/raw.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out $output/raw.QC

# Calculating heterozygosity
plink2 \
    --bfile $output/raw \
    --memory 16000 \
 	  --threads 8 \
    --extract $output/raw.QC.prune.in \
    --keep raw.QC.fam \
    --het \
    --out $output/raw.QC

# Remove '#' from the het file
sed -i "s/#//1" $output/raw.QC.het

# Removing samples with very high or low heterozygosity and keeping valid samples
# Rscript $scripts/sample_het.r $output/

# Alternative SampleQC
plink2 --bfile $output/raw --extract $output/raw.QC.snplist --keep-fam ../../../imputation/usedinpca.txt --write-samples --out $output/sampleQC

# Calling R script to deal with mismatching SNPs
Rscript $scripts/mismatch.r

# Generating final QC'ed target data file
plink2 \
    --bfile $output/raw \
    --make-bed \
    --memory 15000 \
    --threads 8 \
    --keep $output/sampleQC.id \
    --out $output/raw.QC \
    --extract $output/raw.QC.snplist \
    --exclude $output/raw.mismatch \
    --a1-allele $output/raw.a1

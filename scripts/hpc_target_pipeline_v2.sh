#!/bin/bash

helpFunction()
{
	echo ""
	echo "Usage: $0 -i input -b base -p phenotype -o output"
	echo -e "\t-i Path to imputation files"
	echo -e "\t-b Path to QC'd base dataset files"
	echo -e "\t-p Phenotype abbreviation(e.g., SCZ for Schizophrenia)"
	echo -e "\t-o Output directory folder"
	exit 1 #Exit script after printing help
}

while getopts "i:b:p:o:" opt
do
	case "$opt" in
		i ) imp_file_dir="$OPTARG" ;;
		b ) base="$OPTARG" ;;
		p ) phenotype="$OPTARG" ;;
		o ) output="$OPTARG" ;;
		? ) helpFunction ;; #Print helpFunction in case parameter is non-existent
	esac
done

# Print helpFunction in case parameters are emplty
if [ -z "$imp_file_dir" ] || [ -z "$base" ] || [ -z "$phenotype" ] || [ -z "$output" ]
then
	echo "Some or all of the parameters are empty"
	helpFunction
fi

# Extracting the base SNPs from target dataset
for i in {20..22}
do
	bgenix -g "$imp_file_dir"ukb22828_c${i}_b0_v3.bgen -incl-rsids "$base""$phenotype"_rsid.txt > "$output"chr${i}.bgen 
done

#Concatenating the extreacted bgen files
cat-bgen -g "$output"chr*.bgen -og "$output"initial_chr.bgen -clobber

# Writing the index
bgenix -g "$output"initial_chr.bgen -index -clobber

# Removing individual bgen files
for i in {20..22}
do
	rm "$output"chr${i}.bgen
done

# Import the betas into the sqlite database as a table called Betas
sqlite3 "$output"initial_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," "$output"initial_chr.bgen.bgi ".import "$base""$phenotype"_QC.csv Betas"

sqlite3 "$output"/initial_chr.bgen.bgi "DROP TABLE IF EXISTS Joined"
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
 SELECT Variant.*, Betas.CHR, -Betas.BETA FROM Variant INNER JOIN Betas
   ON Variant.chromosome = printf('%02d', Betas.CHR)
   AND Variant.position = Betas.POS
   AND Variant.allele1 = Betas.A1 AND
   Variant.allele2 = Betas.A2;"

# Filter the .bgen file to include only the alleles specified in the Betas for each SNP
bgenix -g "$output"initial_chr.bgen -table Joined  > "$output"single_allelic.bgen

# And produce an index file for the new .bgen
bgenix -g "$output"single_allelic.bgen -index

# Converting file to bed, bim and fam format
plink2 \
	--bgen "$output"single_allelic.bgen ref-first \
 	--sample "$imp_file_dir"ukb22828_c22_b0_v3.sample \
 	--memory 16000 \
	--freq \
 	--make-bed \
 	--out "$output"raw

# Marking ambigous SNPs from target data (I don't think step is relevant)
awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' "$output"raw.afreq > "$output"exclrsIDs_ambiguous.txt

# Standard GWAS QC
plink2 --bfile "$output"raw \
  --memory 16000 \
	--exclude "$output"exclrsIDs_ambiguous.txt \
	--extract-col-cond "$imp_file_dir"ukb_mfi_all_v3.nodup.tsv 9 3 --extract-col-cond-min 0.4 \
	--maf 0.01 \
	--hwe 1e-6 \
	--geno 0.01 \
	--mind 0.01 \
	--rm-dup "exclude-mismatch" \
	--write-snplist \
	--make-just-fam \
	--out "$output"raw.QC

# Pruning to remove highly correlated SNPs
plink2 \
    --bfile "$output"raw \
    --memory 16000 \
    --keep "$output"raw.QC.fam \
    --extract "$output"raw.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out "$output"raw.QC

# Calculating heterozygosity
plink2 \
    --bfile "$output"raw \
    --memory 16000 \
    --extract "$output"raw.QC.prune.in \
    --keep raw.QC.fam \
    --het \
    --out "$output"raw.QC

# Remove '#' from the het file
sed -i "s/#//1" "$output"raw.QC.het

# Alternative SampleQC
plink2 --bfile "$output"raw \
	--extract "$output"raw.QC.snplist \
	--keep-fam "$imp_file_dir"usedinpca.txt \
	--write-samples \
	--out "$output"sampleQC

# Calling R script to deal with mismatching SNPs
Rscript ./mismatch.r "$base""$phenotype" "$output"

# Generating final QC'ed target data file
plink2 \
    --bfile "$output"raw \
    --make-bed \
    --memory 16000 \
    --keep "$output"sampleQC.id \
    --out "$output"raw_QC \
    --extract "$output"raw.QC.snplist \
    --exclude "$output"raw.mismatch \
    --a1-allele "$output"raw.a1
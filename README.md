<code style="color : red">**Not the final README. Created to share and receive feedback.**</code>

## Introduction
The pipeline calculates the Polygenic Score (PGS) using PRS-CS(Ge et al., 2019). It is implemented on bash and is optimized for calculating PGS for UK Biobank data. The functionality is based on guidelines from documentation by [Jennifer Collister and Xiaonan Liu for PGS calculation on UK Biobank data](https://2cjenn.github.io/PRS_Pipeline/) and [Shing Wan Choi's PRS tutorial](https://choishingwan.github.io/PRS-Tutorial/)

## Dependencies
1. [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
2. [R](https://www.r-project.org/)
3. [Plink](https://www.cog-genomics.org/plink/)
4. [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)

*****Python Packages*****

The Python packages required for the PRS-CS calculation are mentioned in [requirements.txt](https://github.com/DuttaPriyatam/PRScs_pipe/blob/master/requirements.txt)

*****R packages*****
1. data.table==1.15.2
2. magrittr==2.0.3

## Getting Started
- Download and install [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and [R](https://www.r-project.org/)
- Download the necessary Python packages
```bash
pip install -r requirements.txt
```
- Download or clone the PRScs_pipe repository
```bash
git clone https://github.com/DuttaPriyatam/PRScs_pipe.git
````
- Create the following directories
1. data - will contain all the relevant input and output data files including base, target, and reference/imputation/sample data
2. apps - to clone the PRS-CS repository
```bash
mkdir data apps
```
- Download or clone the [PRS-CS](https://github.com/getian107/PRScs) git repository inside the apps folder
```bash
git clone https://github.com/getian107/PRScs.git
```
- Build the folder structure for the data folder. From the root directory ./ run the commands mentioned below. Replace "name_of_your_phenotype" with your phenotype of interest. E.g., SCZ for Schizophrenia. 
```bash
mkdir data/imputation data/reference
mkdir data/reference/LD_reference data/reference/ukbiobank_sample
mkdir data/"name_of_your_phenotype"
mkdir data/"name_of_your_phenotype"/base data/"name_of_your_phenotype"/target data/"name_of_your_phenotype"/results
```
- Download the LD reference panel from https://github.com/getian107/PRScs and extract it on the ./data/LD_reference/ directory
- Download your base summary statistics and place it on the ./data/"name_of_your_phenotype"/base/ directory
- Download [gfetch](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=668) and store it in apps directory.
- Go to /data/imputation/ folder and run the following commands to download imputation data(~2.1 TB data).

<code style="color : red">Note - Copy the ukbkey to the working directory to use gfetch in that directory. For further information on ukbkey and gfetch usage please refer to https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Access_Guide.pdf</code>
```bash
../../apps/gfetch 22828 -c1-22
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_bgi.tgz
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_mfi.tgz
```
Extract the .tgz files. 

<code style="color : red">Note - Delete sex chromosomes if not required.</code>

- Download a sample UK Biobank file in the /data/reference/uk_biobank_sample/ folder.
```bash
../../apps/gfetch 22828 -c22 -m
```

## Usage
The pipeline is implemented in 3 steps/scripts:-
1. SNP_QC.sh - Facilitates the preparation and QC of the base summary stats file for PGS calculation
2. hpc_target_pipeline.sh - Works on the UK Biobank imputation data, extracts the SNPs, and does QC of the target dataset
3. calculate_prs.sh - Actually calculates the PGS

**I. Preparing Base Data**

Command-
```bash
SNP_QC.sh -i /path/to/base_summary_statistics -o /path/to/output
```
-i Input file with the path to the directory

-o Output file with phenotype Abbreviation(e.g., /path/to/dir/SCZ for Schizophrenia)

The script will generate the following files:-

i. "name_of_your_phentotyep"_QC.vcf.tsv - QC'd base summary statistics file 

ii. "name_of_your_phentotyep"_for_PRScs.vcf.tsv - Formatted file containing the columns needed for PRScs calculation

iii. "name_of_your_phentotyep"_rsid.txt - rsid list after QC. This will be necessary to extract SNPs from the target imputation data

**II. Preparing Target Data**

Command-
```bash
bash hpc_target_pipeline_v2.sh -i /path/to/impuataion_files -b path/to/QCd_base_summ_stats -p phenotype_abbreviation -o path/to/output_dir
```
-i Path to imputation files

-b Path to QC'd base dataset files

-p Phenotype abbreviation(e.g., SCZ for Schizophrenia)

-o Output directory folder

The script will generate the following files:-

"name_of_your_phenotype"_target_QC.bed, .bim and .fam files - These files will be used for PGS calculation

**III. PGS Calculation**

The script utilizes the PRS-CS repository and plink to generate the PGS. The script generates the PGS for an array of phi values - auto, 1, 1e-2, 1e-4, 1e-6

Command-
```bash
bash calculate_prs.sh -b /path/to/base_PRScs_file -t path/to/target -r path/to/LD_reference -o path/to/output_dir
```
-b Path to the PRS-CS prepared base summary statistics file

-t Path to target bim file(only provide till the prefix for the QCd .bim file) 

-r Path to the LD Reference panel

-p Phenotype abbreviation(e.g., SCZ for Schizophrenia)

-n Total GWAS participants

-o Output directory folder

The script will generate the .score files for each phi parameter used to compute the PGS.

<code style="color : red">**At this point only using auto to calculate the score**</code>


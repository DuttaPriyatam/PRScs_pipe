#!/usr/bin/env Rscript

# Packages to load
packages = c("data.table", "magrittr")

#Checking if package is present and installing/loading the packages
package.check <- lapply(
    packages,
    FUN = function(x) {
        if(!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)

#base <- "/home/pdutta/scratch/data/SCZ/base/"
#target <- "/home/pdutta/scratch/data/SCZ/target/filter_1/"
args <- commandArgs(trailingOnly = TRUE)
base <- args[1]
target <- args[2]

# magrittr allow us to do piping, which help to reduce the amount of intermediate data types
# library(data.table)
# library(magrittr)
# Read in bim file 
bim <- fread(paste(target,"raw.bim", sep="")) %>%
    # Note: . represents the output from previous step
    # The syntax here means, setnames of the data read from
    # the bim file, and replace the original column names by 
    # the new names
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    # And immediately change the alleles to upper cases
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in summary statistic data (require data.table v1.12.0+)
trait <- fread(paste(base,"_QC.vcf.tsv", sep="")) %>%
    # And immediately change the alleles to upper cases
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
# Read in QCed SNPs
qc <- fread(paste(target, "raw.QC.snplist", sep=""), header=F)

# Merge summary statistic with target
info <- merge(bim, trait, by=c("SNP", "CHR", "BP")) %>%
    # And filter out QCed SNPs
    .[SNP %in% qc[,V1]]

# Function for calculating the complementary allele
complement <- function(x){
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 &
                    sapply(B.A2, complement) == A2, SNP]
# Now update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
            sapply(B.A2, complement))]

# identify SNPs that need recoding
recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
# Update the bim file
bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
        list(B.A2, B.A1)]

# identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 &
                    sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
            sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], paste(target, "raw.a1", sep=""), col.names=F, sep="\t")

mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
write.table(mismatch, paste(target, "raw.mismatch", sep=""), quote=F, row.names=F, col.names=F)

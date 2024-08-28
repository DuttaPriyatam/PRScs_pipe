# Taking the input file as as arguement
args <- commandArgs(trailingOnly = TRUE)
base <- args[1]

# Reading the QC file onto R
data <- read.table(file= base, header = TRUE, sep="\t")

# Finding the OR vaues for beta and storing it in the same column 
data$BETA <- exp(data$BETA)

# Renaming the beta column to OR column
names(data)[names(data) == "BETA"] <- "OR"

# Correcting the Standard errors with respect to OR
data$SE <- round(data$OR*data$SE, digit=4)

# Rewriting the QCd base summ stat file
write.table(data, base, quote=F, row.names=F, sep="\t")

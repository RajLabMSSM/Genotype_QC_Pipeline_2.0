# make ancestry PC file for tensorQTL
# using somalier ancestry outputs

# Jack Humphrey 2021


# inputs
# somalier ancestry file

# outputs
# wide-format matrix of the 5 principal components in each donor
library(optparse)
option_list <- list(
    make_option(c('--inFile'), help = "somalier ancestry file", default = ""),
    make_option(c('--outFile'), help = "wide-format PC matrix file", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


inFile <- opt$inFile
outFile <- opt$outFile


# TESTING
#inFile <- "/sc/arion/projects/als-omics/microglia_isoseq/Genotype_QC_Pipeline_2.0/output/roussos/population/roussos_microglia.somalier-ancestry.somalier-ancestry.tsv" 

library(tidyverse)

df <- read_tsv(inFile)

# remove 1000G samples
df <- filter(df, is.na(given_ancestry) )

names(df)[1] <- "donor"

df <- select(df, donor, starts_with("PC") )

# transpose
df_t <- df %>%
    distinct() %>%
    pivot_longer(names_to = "ID", values_to = "value", cols = !donor) %>%
    pivot_wider(names_from = "donor", values_from = "value")

write_tsv(df_t, outFile)



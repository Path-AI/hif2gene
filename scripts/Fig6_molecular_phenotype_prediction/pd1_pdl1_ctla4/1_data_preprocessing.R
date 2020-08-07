library(devtools)
library(readxl)  
library(data.table)
library(stringr)
library(ggplot2)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Source previous custom functions
source("scripts/helper_functions/read_excel.R")

# Load in HIF dataframes
brca <- read_excel_allsheets("data/hifs/brca_hifs.xlsx")
skcm <- read_excel_allsheets("data/hifs/skcm_hifs.xlsx")
stad <- read_excel_allsheets("data/hifs/stad_hifs.xlsx")
lusc <- read_excel_allsheets("data/hifs/lusc_hifs.xlsx")
luad <- read_excel_allsheets("data/hifs/luad_hifs.xlsx")

# brca.defs <- brca$`COLUMN INFO` # feature information
# Extract HIF excel sheet
brca.hifs <- brca$FEATURE 
skcm.hifs <-  skcm$FEATURE
stad.hifs <- stad$FEATURE
lusc.hifs <- lusc$FEATURE
luad.hifs <- luad$FEATURE

dim(brca.hifs)
dim(skcm.hifs)
dim(stad.hifs)
dim(lusc.hifs)
dim(luad.hifs)

sample_count <- dim(brca.hifs)[1] + dim(skcm.hifs)[1] + dim(stad.hifs)[1] + dim(lusc.hifs)[1] + dim(luad.hifs)[1]

# Load in PanImmune dataset
pan.immune <- read.table("data/genetic_data/PanImmune/source/Scores_160_Signatures.tsv", sep="\t", header=TRUE) # full dataset

# Remove metadata columns (currently not relevant)
pan.immune <- pan.immune[,-c(1:2)]

# Transpose PanImmune dataset rows to columns
pan.immune.t <- as.data.frame(t(pan.immune))
pan.immune.t <- setDT(pan.immune.t, keep.rownames = TRUE)
p <- ggplot(pan.immune.t, aes(x=PD1_data)) + geom_histogram(color='black', fill='darkgrey', bins=50)
p + geom_vline(aes(xintercept=0.0), color="red", linetype="dashed", size=1)

# Convert TCGA barcode to interoperable 4-2-4 digit format
pan.immune.t <- pan.immune.t %>% dplyr::rename(bcr_patient_barcode = rn)
pan.immune.t <- dplyr::mutate(
  pan.immune.t, bcr_patient_barcode = substring(bcr_patient_barcode, 1, 12))
pan.immune.t <- dplyr::mutate(
  pan.immune.t, bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "[.]", "-"))

# Filter PanImmune dataset for each cancer type
# BRCA
brca.barcodes <- brca.hifs$case
pan.immune.brca <- dplyr::filter(pan.immune.t, bcr_patient_barcode %in% brca.barcodes)
brca.barcodes[duplicated(brca.barcodes)] # Duplicated barcodes: same patient, different samples
cat(sprintf("Pct Match (Patient-Level) HIF & PanImmune=%i/%i=%f\n", 
            n_distinct(pan.immune.brca$bcr_patient_barcode), 
            n_distinct(brca.barcodes),
            n_distinct(pan.immune.brca$bcr_patient_barcode) / n_distinct(brca.barcodes) # Percent matched between HIF and MAF datasets (Patient-Level)
  )
)
saveRDS(pan.immune.brca, file="data/genetic_data/PanImmune/processed/pan.immune.brca.rds") # File containing only BRCA barcode-matched data
temp <- distinct(pan.immune.brca, bcr_patient_barcode)
brca.hifs.with.panimmune.data <- dplyr::filter(brca.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & MAF=%i/%i=%f\n", 
            dim(brca.hifs.with.panimmune.data)[1], 
            dim(brca.hifs)[1],
            dim(brca.hifs.with.panimmune.data)[1] / dim(brca.hifs)[1] # Percent matched between HIF and PanImmune datasets (Sample-Level)
  )
)
# Patient barcodes not found in the PanImmune dataset are "missing" (unknown values) --> filter out
saveRDS(brca.hifs.with.panimmune.data, file="data/genetic_data/PanImmune/processed/brca.hifs.with.panimmune.data.rds") # BRCA HIF dataset filtered for patients / samples with non-missing PanImmune data present

# SKCM
skcm.barcodes <- skcm.hifs$case
pan.immune.skcm <- dplyr::filter(pan.immune.t, bcr_patient_barcode %in% skcm.barcodes)
skcm.barcodes[duplicated(skcm.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & PanImmune=%i/%i=%f\n", 
            n_distinct(pan.immune.skcm$bcr_patient_barcode), 
            n_distinct(skcm.barcodes),
            n_distinct(pan.immune.skcm$bcr_patient_barcode) / n_distinct(skcm.barcodes)
  )
)
saveRDS(pan.immune.skcm, file="data/genetic_data/PanImmune/processed/pan.immune.skcm.rds")
temp <- distinct(pan.immune.skcm, bcr_patient_barcode)
skcm.hifs.with.panimmune.data <- dplyr::filter(skcm.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & MAF=%i/%i=%f\n", 
            dim(skcm.hifs.with.panimmune.data)[1], 
            dim(skcm.hifs)[1],
            dim(skcm.hifs.with.panimmune.data)[1] / dim(skcm.hifs)[1] 
  )
)
saveRDS(skcm.hifs.with.panimmune.data, file="data/genetic_data/PanImmune/processed/skcm.hifs.with.panimmune.data.rds") 

# STAD
stad.barcodes <- stad.hifs$case
pan.immune.stad <- dplyr::filter(pan.immune.t, bcr_patient_barcode %in% stad.barcodes)
stad.barcodes[duplicated(stad.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & PanImmune=%i/%i=%f\n", 
            n_distinct(pan.immune.stad$bcr_patient_barcode), 
            n_distinct(stad.barcodes),
            n_distinct(pan.immune.stad$bcr_patient_barcode) / n_distinct(stad.barcodes)
  )
)
saveRDS(pan.immune.stad, file="data/genetic_data/PanImmune/processed/pan.immune.stad.rds")
temp <- distinct(pan.immune.stad, bcr_patient_barcode)
stad.hifs.with.panimmune.data <- dplyr::filter(stad.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & MAF=%i/%i=%f\n", 
            dim(stad.hifs.with.panimmune.data)[1], 
            dim(stad.hifs)[1],
            dim(stad.hifs.with.panimmune.data)[1] / dim(stad.hifs)[1] 
  )
)
saveRDS(stad.hifs.with.panimmune.data, file="data/genetic_data/PanImmune/processed/stad.hifs.with.panimmune.data.rds") 

# LUSC
lusc.barcodes <- lusc.hifs$case
pan.immune.lusc <- dplyr::filter(pan.immune.t, bcr_patient_barcode %in% lusc.barcodes)
lusc.barcodes[duplicated(lusc.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & PanImmune=%i/%i=%f\n", 
            n_distinct(pan.immune.lusc$bcr_patient_barcode), 
            n_distinct(lusc.barcodes),
            n_distinct(pan.immune.lusc$bcr_patient_barcode) / n_distinct(lusc.barcodes)
  )
)
saveRDS(pan.immune.lusc, file="data/genetic_data/PanImmune/processed/pan.immune.lusc.rds")
temp <- distinct(pan.immune.lusc, bcr_patient_barcode)
lusc.hifs.with.panimmune.data <- dplyr::filter(lusc.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & MAF=%i/%i=%f\n", 
            dim(lusc.hifs.with.panimmune.data)[1], 
            dim(lusc.hifs)[1],
            dim(lusc.hifs.with.panimmune.data)[1] / dim(lusc.hifs)[1] 
  )
)
saveRDS(lusc.hifs.with.panimmune.data, file="data/genetic_data/PanImmune/processed/lusc.hifs.with.panimmune.data.rds") 

# LUAD
luad.barcodes <- luad.hifs$case
pan.immune.luad <- dplyr::filter(pan.immune.t, bcr_patient_barcode %in% luad.barcodes)
luad.barcodes[duplicated(luad.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & PanImmune=%i/%i=%f\n", 
            n_distinct(pan.immune.luad$bcr_patient_barcode), 
            n_distinct(luad.barcodes),
            n_distinct(pan.immune.luad$bcr_patient_barcode) / n_distinct(luad.barcodes)
  )
)
saveRDS(pan.immune.luad, file="data/genetic_data/PanImmune/processed/pan.immune.luad.rds")
temp <- distinct(pan.immune.luad, bcr_patient_barcode)
luad.hifs.with.panimmune.data <- dplyr::filter(luad.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & MAF=%i/%i=%f\n", 
            dim(luad.hifs.with.panimmune.data)[1], 
            dim(luad.hifs)[1],
            dim(luad.hifs.with.panimmune.data)[1] / dim(luad.hifs)[1] 
  )
)
saveRDS(luad.hifs.with.panimmune.data, file="data/genetic_data/PanImmune/processed/luad.hifs.with.panimmune.data.rds") 

# Plot histograms
nbins <- 150
p <- ggplot() + geom_histogram(aes(x=pan.immune.t$PD1_data), color='black', fill='darkgrey', bins=nbins)
p <- p + geom_histogram(aes(x=pan.immune.brca$PD1_data), color='black', fill='red', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.luad$PD1_data), color='black', fill='blue', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.lusc$PD1_data), color='black', fill='green', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.stad$PD1_data), color='black', fill='yellow', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.skcm$PD1_data), color='black', fill='orange', bins=nbins, alpha=0.5)
p + xlab("PD-1 Score")

nbins <- 150
p <- ggplot() + geom_histogram(aes(x=pan.immune.t$PDL1_data), color='black', fill='darkgrey', bins=nbins)
p <- p + geom_histogram(aes(x=pan.immune.brca$PDL1_data), color='black', fill='red', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.luad$PDL1_data), color='black', fill='blue', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.lusc$PDL1_data), color='black', fill='green', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.stad$PDL1_data), color='black', fill='yellow', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.skcm$PDL1_data), color='black', fill='orange', bins=nbins, alpha=0.5)
p + xlab("PDL-1 Score")

nbins <- 150
p <- ggplot() + geom_histogram(aes(x=pan.immune.t$CTLA4_data), color='black', fill='darkgrey', bins=nbins)
p <- p + geom_histogram(aes(x=pan.immune.brca$CTLA4_data), color='black', fill='red', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.luad$CTLA4_data), color='black', fill='blue', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.lusc$CTLA4_data), color='black', fill='green', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.stad$CTLA4_data), color='black', fill='yellow', bins=nbins, alpha=0.5)
p <- p + geom_histogram(aes(x=pan.immune.skcm$CTLA4_data), color='black', fill='orange', bins=nbins, alpha=0.5)
p + xlab("CTLA-4 Score")


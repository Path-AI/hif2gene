library(devtools)
library(readxl)  
library(data.table)
library(stringr)

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

# Load in DDR scores for TCGA data
ddr.scores.full <- read_excel_allsheets("data/genetic_data/HRD/source/TCGA_DDR_Data_Resources.xlsx") # full dataset
ddr.scores <- ddr.scores.full$`DDR footprints`

# Remove metadata rows
names(ddr.scores) <- ddr.scores[3,]
ddr.scores <- ddr.scores[-c(1,2,3),]

# Remove metadata columns (currently not relevant)
ddr.scores <- ddr.scores[,-c(2,3,4)]
ddr.scores <- ddr.scores %>% dplyr::rename(bcr_patient_barcode = patient_barcode)

# Filter PanImmune dataset for each cancer type
# BRCA
brca.barcodes <- brca.hifs$case
ddr.scores.brca <- dplyr::filter(ddr.scores, bcr_patient_barcode %in% brca.barcodes)
brca.barcodes[duplicated(brca.barcodes)] # Duplicated barcodes: same patient, different samples
cat(sprintf("Pct Match (Patient-Level) HIF & DDR Scores=%i/%i=%f\n", 
            n_distinct(ddr.scores.brca$bcr_patient_barcode), 
            n_distinct(brca.barcodes),
            n_distinct(ddr.scores.brca$bcr_patient_barcode) / n_distinct(brca.barcodes) # Percent matched between HIF and DDR scores (Patient-Level)
  )
)
saveRDS(ddr.scores.brca, file="data/genetic_data/HRD/processed/ddr.scores.brca.rds") # File containing only BRCA barcode-matched data
temp <- distinct(ddr.scores.brca, bcr_patient_barcode)
brca.hifs.with.ddr.scores <- dplyr::filter(brca.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & DDR Scores=%i/%i=%f\n", 
            dim(brca.hifs.with.ddr.scores)[1], 
            dim(brca.hifs)[1],
            dim(brca.hifs.with.ddr.scores)[1] / dim(brca.hifs)[1] # Percent matched between HIF and DDR scores (Sample-Level)
  )
)
# Patient barcodes not found in the PanImmune dataset are "missing" (unknown values) --> filter out
saveRDS(brca.hifs.with.ddr.scores, file="data/genetic_data/HRD/processed/brca.hifs.with.ddr.scores.rds") # BRCA HIF dataset filtered for patients / samples with non-missing TCGA DDR data present

# SKCM
skcm.barcodes <- skcm.hifs$case
ddr.scores.skcm <- dplyr::filter(ddr.scores, bcr_patient_barcode %in% skcm.barcodes)
skcm.barcodes[duplicated(skcm.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & DDR Scores=%i/%i=%f\n", 
            n_distinct(ddr.scores.skcm$bcr_patient_barcode), 
            n_distinct(skcm.barcodes),
            n_distinct(ddr.scores.skcm$bcr_patient_barcode) / n_distinct(skcm.barcodes) 
  )
)
saveRDS(ddr.scores.skcm, file="data/genetic_data/HRD/processed/skcm/ddr.scores.skcm.rds")
temp <- distinct(ddr.scores.skcm, bcr_patient_barcode)
skcm.hifs.with.ddr.scores <- dplyr::filter(skcm.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & DDR Scores=%i/%i=%f\n", 
            dim(skcm.hifs.with.ddr.scores)[1], 
            dim(skcm.hifs)[1],
            dim(skcm.hifs.with.ddr.scores)[1] / dim(skcm.hifs)[1] 
  )
)
saveRDS(skcm.hifs.with.ddr.scores, file="data/genetic_data/HRD/processed/skcm.hifs.with.ddr.scores.rds") 

# STAD
stad.barcodes <- stad.hifs$case
ddr.scores.stad <- dplyr::filter(ddr.scores, bcr_patient_barcode %in% stad.barcodes)
stad.barcodes[duplicated(stad.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & DDR Scores=%i/%i=%f\n", 
            n_distinct(ddr.scores.stad$bcr_patient_barcode), 
            n_distinct(stad.barcodes),
            n_distinct(ddr.scores.stad$bcr_patient_barcode) / n_distinct(stad.barcodes) 
  )
)
saveRDS(ddr.scores.stad, file="data/genetic_data/HRD/processed/ddr.scores.stad.rds")
temp <- distinct(ddr.scores.stad, bcr_patient_barcode)
stad.hifs.with.ddr.scores <- dplyr::filter(stad.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & DDR Scores=%i/%i=%f\n", 
            dim(stad.hifs.with.ddr.scores)[1], 
            dim(stad.hifs)[1],
            dim(stad.hifs.with.ddr.scores)[1] / dim(stad.hifs)[1] 
  )
)
saveRDS(stad.hifs.with.ddr.scores, file="data/genetic_data/HRD/processed/stad.hifs.with.ddr.scores.rds") 

# LUSC
lusc.barcodes <- lusc.hifs$case
ddr.scores.lusc <- dplyr::filter(ddr.scores, bcr_patient_barcode %in% lusc.barcodes)
lusc.barcodes[duplicated(lusc.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & DDR Scores=%i/%i=%f\n", 
            n_distinct(ddr.scores.lusc$bcr_patient_barcode), 
            n_distinct(lusc.barcodes),
            n_distinct(ddr.scores.lusc$bcr_patient_barcode) / n_distinct(lusc.barcodes) 
  )
)
saveRDS(ddr.scores.lusc, file="data/genetic_data/HRD/processed/ddr.scores.lusc.rds")
temp <- distinct(ddr.scores.lusc, bcr_patient_barcode)
lusc.hifs.with.ddr.scores <- dplyr::filter(lusc.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & DDR Scores=%i/%i=%f\n", 
            dim(lusc.hifs.with.ddr.scores)[1], 
            dim(lusc.hifs)[1],
            dim(lusc.hifs.with.ddr.scores)[1] / dim(lusc.hifs)[1] 
  )
)
saveRDS(lusc.hifs.with.ddr.scores, file="data/genetic_data/HRD/processed/lusc.hifs.with.ddr.scores.rds") 

# LUAD
luad.barcodes <- luad.hifs$case
ddr.scores.luad <- dplyr::filter(ddr.scores, bcr_patient_barcode %in% luad.barcodes)
luad.barcodes[duplicated(luad.barcodes)] 
cat(sprintf("Pct Match (Patient-Level) HIF & DDR Scores=%i/%i=%f\n", 
            n_distinct(ddr.scores.luad$bcr_patient_barcode), 
            n_distinct(luad.barcodes),
            n_distinct(ddr.scores.luad$bcr_patient_barcode) / n_distinct(luad.barcodes) 
  )
)
saveRDS(ddr.scores.luad, file="data/genetic_data/HRD/processed/ddr.scores.luad.rds")
temp <- distinct(ddr.scores.luad, bcr_patient_barcode)
luad.hifs.with.ddr.scores <- dplyr::filter(luad.hifs, case %in% temp$bcr_patient_barcode)
cat(sprintf("Pct Match (Sample-Level) HIF & DDR Scores=%i/%i=%f\n", 
            dim(luad.hifs.with.ddr.scores)[1], 
            dim(luad.hifs)[1],
            dim(luad.hifs.with.ddr.scores)[1] / dim(luad.hifs)[1] 
  )
)
saveRDS(luad.hifs.with.ddr.scores, file="data/genetic_data/HRD/processed/luad.hifs.with.ddr.scores.rds") 

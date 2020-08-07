library(devtools)
library(readxl)  

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Load in HIF dataframes (filtered for samples / patients with genotype data)
brca.hifs.with.panimmune.data <- readRDS(file="data/genetic_data/PanImmune/processed/brca.hifs.with.panimmune.data.rds")
brca.hifs.with.panimmune.data <- dplyr::select(brca.hifs.with.panimmune.data, bcr_patient_barcode, everything()) # Key of interest = TCGA barcode, move to first column of dataframe

# Load in filtered (by TCGA barcodes in PanImmune dataset) MAF dataframes 
pan.immune.brca <- readRDS(file="data/genetic_data/PanImmune/processed/pan.immune.brca.rds")

# (Optional) 
# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
brca.barcodes <- brca.hifs.with.panimmune.data$bcr_patient_barcode 
length(brca.barcodes[duplicated(brca.barcodes)]) # Duplicated barcodes: same patient, different samples
brca.hifs.with.panimmune.data <- brca.hifs.with.panimmune.data %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(brca.hifs.with.panimmune.data)
# Second filter, if tumor areas are equal, choose randomly
brca.hifs.with.panimmune.data <- brca.hifs.with.panimmune.data %>% group_by(bcr_patient_barcode) %>% sample_n(1)
dim(brca.hifs.with.panimmune.data)

# Join HIF and PanImmune datasets
joined <- dplyr::left_join(brca.hifs.with.panimmune.data, pan.immune.brca, by="bcr_patient_barcode")

# Saved joined dataframe
saveRDS(joined, file="data/genetic_data/PanImmune/processed/brca.pan.immune.joined.rds")
write.csv(joined, file="data/genetic_data/PanImmune/processed/brca.pan.immune.joined.csv")

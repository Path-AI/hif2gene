library(devtools)
library(readxl)  
library(caret)
library(mclust)
library(dplyr)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Load in HIF dataframes (filtered for samples / patients with genotype data)
brca.hifs.with.tigit.data <- readRDS(file="data/genetic_data/TIGIT/processed/brca.hifs.with.tigit.data.rds")
brca.hifs.with.tigit.data <- dplyr::select(brca.hifs.with.tigit.data, bcr_patient_barcode, everything()) # Key of interest = TCGA barcode, move to first column of dataframe

# Load in filtered (by TCGA barcodes in PanImmune dataset) MAF dataframes 
tigit.brca <- readRDS(file="data/genetic_data/TIGIT/processed/tigit.brca.rds")

# (Optional) 
# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
brca.barcodes <- brca.hifs.with.tigit.data$bcr_patient_barcode 
length(brca.barcodes[duplicated(brca.barcodes)]) # Duplicated barcodes: same patient, different samples
brca.hifs.with.tigit.data <- brca.hifs.with.tigit.data %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(brca.hifs.with.tigit.data)
# Second filter, if tumor areas are equal, choose randomly
brca.hifs.with.tigit.data <- brca.hifs.with.tigit.data %>% group_by(bcr_patient_barcode) %>% sample_n(1)
dim(brca.hifs.with.tigit.data)

# Check dupliates in TIGIT data
tigit.barcodes <- tigit.brca$bcr_patient_barcode
length(tigit.barcodes[duplicated(tigit.barcodes)])

###############################
# Find TIGIT Binary Threshold #
###############################
# Plot TIGIT score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(tigit.brca$TIGIT, G=2, model="V")
summary(gmm.fit, parameters=TRUE)
tigit.brca$GMM_Cluster <- gmm.fit$classification
tigit.brca$GMM_Cluster <- mapvalues(tigit.brca$GMM_Cluster,
                                from=c(1,2,3,4),
                                to=c('r', 'o', 'y', 'b'))

# Manually find intersection between two density plots
tigit.threshold <- 0.25
p <- ggplot(tigit.brca, aes(x=TIGIT)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=TIGIT), alpha=0) + geom_density(aes(x=TIGIT, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=tigit.threshold), color="red", linetype="dashed", size=1)

# After plotting, remove GMM_Cluster column
tigit.brca = tigit.brca[,-length(tigit.brca)[1]]

# Add binary column to TIGIT data
tigit.brca$TIGIT_label <- ifelse(tigit.brca$TIGIT > tigit.threshold, 1, 0)
tigit.brca <- tigit.brca[,c(-2)]

######################

# Reduce TIGIT data to patient-level
# If binary labels disagree, remove both duplicated samples
# If binary labels agree, remove duplicated but keep one
dim(tigit.brca)
tigit.brca.cleaned <- tigit.brca %>% dplyr::group_by(bcr_patient_barcode) %>% dplyr::summarise(TIGIT_label = mean(TIGIT_label))
tigit.brca.cleaned <- tigit.brca.cleaned[tigit.brca.cleaned$TIGIT_label == 1.0 | tigit.brca.cleaned$TIGIT_label == 0.0, ] # labels disagree
dim(tigit.brca.cleaned)
sum(tigit.brca.cleaned$TIGIT_label)/length(tigit.brca.cleaned$TIGIT_label) # Percent above threshold cutoff

# Join HIF and PanImmune datasets
joined <- dplyr::left_join(brca.hifs.with.tigit.data, tigit.brca.cleaned, by="bcr_patient_barcode")

# Extract Tissue Source Site Info (TSS)
joined$TSS <- substr(joined$bcr_patient_barcode, 6, 7)
site.breakdown <- as.data.frame(table(joined$TSS))
colnames(site.breakdown) <- c('TSS', 'Sample Count')

# Extract relevant columns (607 HIFs + 6 genotypes = 613 columns)
joined.relevant <- data.frame(joined[,c(4:610)], 
                              joined$TIGIT_label, joined$TSS)
dim(joined.relevant)

# Throw out rows with high missigness, otherwise mean imputation
na.rows <- joined.relevant[rowSums(is.na(joined.relevant)) > 0,] # two rows comprising all NA values
joined.cleaned <- joined.relevant[-as.numeric(rownames(na.rows)),]
dim(joined.cleaned)
sum(is.na(joined.cleaned)) # verify 0 missingness

# Hold-out TSS's
tss.holdout <- c('BH', 'A2')
use <- dplyr::filter(joined.cleaned, !joined.TSS %in% tss.holdout)
holdout <- dplyr::filter(joined.cleaned, joined.TSS %in% tss.holdout)
dim(use)[1]/dim(joined.cleaned)[1]
dim(holdout)[1]/dim(joined.cleaned)[1]

sum(joined.cleaned$joined.TIGIT_label)/dim(joined.cleaned)[1]
sum(use$joined.TIGIT_label)/dim(use)[1]
sum(holdout$joined.TIGIT_label)/dim(holdout)[1]

# Save split dataframes
tigit.use <- use[,-c(609)]
tigit.hold.out <- holdout[,-c(609)]
write.csv(tigit.use, file="data/genetic_data/TIGIT/datasets/brca.tigit.joined.USE.csv")
saveRDS(tigit.use, file="data/genetic_data/TIGIT/datasets/brca.tigit.joined.USE.rds")
write.csv(tigit.hold.out, file="data/genetic_data/TIGIT/datasets/brca.tigit.joined.HOLDOUT.csv")
saveRDS(tigit.hold.out, file="data/genetic_data/TIGIT/datasets/brca.tigit.joined.HOLDOUT.rds")

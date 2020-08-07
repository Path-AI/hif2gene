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
lusc.hifs.with.tigit.data <- readRDS(file="data/genetic_data/TIGIT/processed/lusc.hifs.with.tigit.data.rds")
lusc.hifs.with.tigit.data <- dplyr::select(lusc.hifs.with.tigit.data, bcr_patient_barcode, everything()) # Key of interest = TCGA barcode, move to first column of dataframe

# Load in filtered (by TCGA barcodes in PanImmune dataset) MAF dataframes 
tigit.lusc <- readRDS(file="data/genetic_data/TIGIT/processed/tigit.lusc.rds")

# (Optional) 
# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
lusc.barcodes <- lusc.hifs.with.tigit.data$bcr_patient_barcode 
length(lusc.barcodes[duplicated(lusc.barcodes)]) # Duplicated barcodes: same patient, different samples
lusc.hifs.with.tigit.data <- lusc.hifs.with.tigit.data %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`))
dim(lusc.hifs.with.tigit.data)
# Second filter, if tumor areas are equal, choose randomly
lusc.hifs.with.tigit.data <- lusc.hifs.with.tigit.data %>% group_by(bcr_patient_barcode) %>% sample_n(1)
dim(lusc.hifs.with.tigit.data)

# Check dupliates in TIGIT data
tigit.barcodes <- tigit.lusc$bcr_patient_barcode
length(tigit.barcodes[duplicated(tigit.barcodes)])

###############################
# Find TIGIT Binary Threshold #
###############################
# Plot TIGIT score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(tigit.lusc$TIGIT, G=4, model="V")
summary(gmm.fit, parameters=TRUE)
tigit.lusc$GMM_Cluster <- gmm.fit$classification
tigit.lusc$GMM_Cluster <- mapvalues(tigit.lusc$GMM_Cluster,
                                from=c(1,2,3,4),
                                to=c('b', 'b', 'b', 'r'))

# Manually find intersection between two density plots
tigit.threshold <- 0.92
p <- ggplot(tigit.lusc, aes(x=TIGIT)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=TIGIT), alpha=0) + geom_density(aes(x=TIGIT, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=tigit.threshold), color="red", linetype="dashed", size=1)

# After plotting, remove GMM_Cluster column
tigit.lusc = tigit.lusc[,-length(tigit.lusc)[1]]

# Add binary column to TIGIT data
tigit.lusc$TIGIT_label <- ifelse(tigit.lusc$TIGIT > tigit.threshold, 1, 0)
tigit.lusc <- tigit.lusc[,c(-2)]

######################

# Reduce TIGIT data to patient-level
# If binary labels disagree, remove both duplicated samples
# If binary labels agree, remove duplicated but keep one
dim(tigit.lusc)
tigit.lusc.cleaned <- tigit.lusc %>% dplyr::group_by(bcr_patient_barcode) %>% dplyr::summarise(TIGIT_label = mean(TIGIT_label))
tigit.lusc.cleaned <- tigit.lusc.cleaned[tigit.lusc.cleaned$TIGIT_label != 0.5, ] # labels disagree
dim(tigit.lusc.cleaned)
sum(tigit.lusc.cleaned$TIGIT_label)/length(tigit.lusc.cleaned$TIGIT_label) # Percent above threshold cutoff

# Join HIF and PanImmune datasets
joined <- dplyr::left_join(lusc.hifs.with.tigit.data, tigit.lusc.cleaned, by="bcr_patient_barcode")

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
tss.holdout <- c('85', '66', '33')
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
write.csv(tigit.use, file="data/genetic_data/TIGIT/datasets/lusc.tigit.joined.USE.csv")
saveRDS(tigit.use, file="data/genetic_data/TIGIT/datasets/lusc.tigit.joined.USE.rds")
write.csv(tigit.hold.out, file="data/genetic_data/TIGIT/datasets/lusc.tigit.joined.HOLDOUT.csv")
saveRDS(tigit.hold.out, file="data/genetic_data/TIGIT/datasets/lusc.tigit.joined.HOLDOUT.rds")

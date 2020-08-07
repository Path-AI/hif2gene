library(devtools)
library(readxl)  
library(caret)
library(mclust)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Load in HIF dataframes (filtered for samples / patients with genotype data)
luad.hifs.with.panimmune.data <- readRDS(file="data/genetic_data/PanImmune/processed/luad.hifs.with.panimmune.data.rds")
luad.hifs.with.panimmune.data <- dplyr::select(luad.hifs.with.panimmune.data, bcr_patient_barcode, everything()) # Key of interest = TCGA barcode, move to first column of dataframe

# Load in filtered (by TCGA barcodes in PanImmune dataset) MAF dataframes 
pan.immune.luad <- readRDS(file="data/genetic_data/PanImmune/processed/pan.immune.luad.rds")

# (Optional) 
# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
luad.barcodes <- luad.hifs.with.panimmune.data$bcr_patient_barcode 
length(luad.barcodes[duplicated(luad.barcodes)]) # Duplicated barcodes: same patient, different samples
luad.hifs.with.panimmune.data <- luad.hifs.with.panimmune.data %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(luad.hifs.with.panimmune.data)
# Second filter, if tumor areas are equal, choose randomly
luad.hifs.with.panimmune.data <- luad.hifs.with.panimmune.data %>% group_by(bcr_patient_barcode) %>% sample_n(1)
dim(luad.hifs.with.panimmune.data)

# Join HIF and PanImmune datasets
joined <- dplyr::left_join(luad.hifs.with.panimmune.data, pan.immune.luad, by="bcr_patient_barcode")

# Outcomes of interest (to predict)
# 'PD1_data', 'PDL1_data', 'CTLA4_data'

#######################
# Find PD-1 Threshold #
#######################
# Plot PD-1 score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(joined$PD1_data, G=2, model="V")
summary(gmm.fit, parameters=TRUE)
joined$GMM_Cluster <- gmm.fit$classification
joined$GMM_Cluster <- mapvalues(joined$GMM_Cluster,
                                        from=c(1,2,3,4),
                                        to=c('r', 'o', 'y', 'b'))

# Manually find intersection between two density plots
# Plot PD-1 score distribution with chosen threshold
PD1.threshold <- 0.87
p <- ggplot(joined, aes(x=PD1_data)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=PD1_data), alpha=0) + geom_density(aes(x=PD1_data, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=PD1.threshold), color="red", linetype="dashed", size=1)
sum(joined$PD1_data > PD1.threshold)/length(joined$PD1_data) # Percent above threshold cutoff

# After plotting, remove GMM_Cluster column
joined = joined[,-length(joined)[1]]

########################
# Find PDL-1 Threshold #
########################
# Plot PDL-1 score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(joined$PDL1_data, G=4, model="V")
summary(gmm.fit, parameters=TRUE)
joined$GMM_Cluster <- gmm.fit$classification
joined$GMM_Cluster <- mapvalues(joined$GMM_Cluster,
                                from=c(1,2,3,4),
                                to=c('r', 'r', 'b', 'b'))

# Manually find intersection between two density plots
# Plot PD-1 score distribution with chosen threshold
PDL1.threshold <- 0.66
p <- ggplot(joined, aes(x=PDL1_data)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=PDL1_data), alpha=0) + geom_density(aes(x=PDL1_data, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=PDL1.threshold), color="red", linetype="dashed", size=1)
sum(joined$PDL1_data > PDL1.threshold)/length(joined$PDL1_data) # Percent above threshold cutoff

# After plotting, remove GMM_Cluster column
joined = joined[,-length(joined)[1]]

#########################
# Find CTLA-4 Threshold #
#########################
# Plot CTLA-4 score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(joined$CTLA4_data, G=4, model="V")
summary(gmm.fit, parameters=TRUE)
joined$GMM_Cluster <- gmm.fit$classification
joined$GMM_Cluster <- mapvalues(joined$GMM_Cluster,
                                from=c(1,2,3,4),
                                to=c('r', 'r', 'b', 'b'))

# Manually find intersection between two density plots
# Plot PD-1 score distribution with chosen threshold
CTLA4.threshold <- 0.97
p <- ggplot(joined, aes(x=CTLA4_data)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=CTLA4_data), alpha=0) + geom_density(aes(x=CTLA4_data, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=CTLA4.threshold), color="red", linetype="dashed", size=1)
sum(joined$CTLA4_data > CTLA4.threshold)/length(joined$CTLA4_data) # Percent above threshold cutoff

# After plotting, remove GMM_Cluster column
joined = joined[,-length(joined)[1]]

######################

# Convert continuous (regression) label --> binary (classification) label
# Above / below Z-score = 0
joined$PD1_data <- ifelse(joined$PD1_data > PD1.threshold, 1, 0)
joined$PDL1_data <- ifelse(joined$PDL1_data > PDL1.threshold, 1, 0)
joined$CTLA4_data <- ifelse(joined$CTLA4_data > CTLA4.threshold, 1, 0)

# Extract Tissue Source Site Info (TSS)
joined$TSS <- substr(joined$bcr_patient_barcode, 6, 7)
site.breakdown <- as.data.frame(table(joined$TSS))
colnames(site.breakdown) <- c('TSS', 'Sample Count')

# Extract relevant columns (607 HIFs + 6 genotypes = 613 columns)
joined.relevant <- data.frame(joined[,c(4:610)], 
                              joined$PD1_data, joined$PDL1_data, joined$CTLA4_data, joined$TSS)
dim(joined.relevant)

# Throw out rows with high missigness, otherwise mean imputation
na.rows <- joined.relevant[rowSums(is.na(joined.relevant)) > 0,] # two rows comprising all NA values
joined.cleaned <- joined.relevant[-as.numeric(rownames(na.rows)),]
dim(joined.cleaned)
sum(is.na(joined.cleaned)) # verify 0 missingness

# Hold-out TSS's
tss.holdout <- c('55', '50')
use <- dplyr::filter(joined.cleaned, !joined.TSS %in% tss.holdout)
holdout <- dplyr::filter(joined.cleaned, joined.TSS %in% tss.holdout)
dim(use)[1]/dim(joined.cleaned)[1]
dim(holdout)[1]/dim(joined.cleaned)[1]

sum(joined.cleaned$joined.PD1_data)/dim(joined.cleaned)[1]
sum(use$joined.PD1_data)/dim(use)[1]
sum(holdout$joined.PD1_data)/dim(holdout)[1]

sum(joined.cleaned$joined.PDL1_data)/dim(joined.cleaned)[1]
sum(use$joined.PDL1_data)/dim(use)[1]
sum(holdout$joined.PDL1_data)/dim(holdout)[1]

sum(joined.cleaned$joined.CTLA4_data)/dim(joined.cleaned)[1]
sum(use$joined.CTLA4_data)/dim(use)[1]
sum(holdout$joined.CTLA4_data)/dim(holdout)[1]

# Save split dataframes
pd1.use <- use[,-c(609, 610, 611)]
pd1.hold.out <- holdout[,-c(609, 610, 611)]
write.csv(pd1.use, file="data/genetic_data/PanImmune/datasets/luad.PD1.joined.USE.csv")
saveRDS(pd1.use, file="data/genetic_data/PanImmune/datasets/luad.PD1.joined.USE.rds")
write.csv(pd1.hold.out, file="data/genetic_data/PanImmune/datasets/luad.PD1.joined.HOLDOUT.csv")
saveRDS(pd1.hold.out, file="data/genetic_data/PanImmune/datasets/luad.PD1.joined.HOLDOUT.rds")

pdl1.use <- use[,-c(608, 610, 611)]
pdl1.hold.out <- holdout[,-c(608, 610, 611)]
write.csv(pdl1.use, file="data/genetic_data/PanImmune/datasets/luad.PDL1.joined.USE.csv")
saveRDS(pdl1.use, file="data/genetic_data/PanImmune/datasets/luad.PDL1.joined.USE.rds")
write.csv(pdl1.hold.out, file="data/genetic_data/PanImmune/datasets/luad.PDL1.joined.HOLDOUT.csv")
saveRDS(pdl1.hold.out, file="data/genetic_data/PanImmune/datasets/luad.PDL1.joined.HOLDOUT.rds")

ctla4.use <- use[,-c(608, 609, 611)]
ctla4.hold.out <- holdout[,-c(608, 609, 611)]
write.csv(ctla4.use, file="data/genetic_data/PanImmune/datasets/luad.CTLA4.joined.USE.csv")
saveRDS(ctla4.use, file="data/genetic_data/PanImmune/datasets/luad.CTLA4.joined.USE.rds")
write.csv(ctla4.hold.out, file="data/genetic_data/PanImmune/datasets/luad.CTLA4.joined.HOLDOUT.csv")
saveRDS(ctla4.hold.out, file="data/genetic_data/PanImmune/datasets/luad.CTLA4.joined.HOLDOUT.rds")

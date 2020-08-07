library(devtools)
library(readxl)  
library(caret)
library(mclust, quietly=TRUE)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Load in HIF dataframes (filtered for samples / patients with TCGA_DDR data)
brca.hifs.with.ddr.scores <- readRDS(file="data/genetic_data/HRD/processed/brca.hifs.with.ddr.scores.rds")
brca.hifs.with.ddr.scores <- dplyr::select(brca.hifs.with.ddr.scores, bcr_patient_barcode, everything()) # Key of interest = TCGA barcode, move to first column of dataframe

# Load in DDR scores data (filtered by TCGA barcodes) 
ddr.scores.brca <- readRDS(file="data/genetic_data/HRD/processed/ddr.scores.brca.rds")

# (Optional) 
# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
brca.barcodes <- brca.hifs.with.ddr.scores$bcr_patient_barcode 
length(brca.barcodes[duplicated(brca.barcodes)]) # Duplicated barcodes: same patient, different samples
brca.hifs.with.ddr.scores <- brca.hifs.with.ddr.scores %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(brca.hifs.with.ddr.scores)
# Second filter, if tumor areas are equal, choose randomly
brca.hifs.with.ddr.scores <- brca.hifs.with.ddr.scores %>% group_by(bcr_patient_barcode) %>% sample_n(1)
dim(brca.hifs.with.ddr.scores)

# Join HIF and DDR datasets
joined <- dplyr::left_join(brca.hifs.with.ddr.scores, ddr.scores.brca, by="bcr_patient_barcode")
joined$HRD_Score <- as.numeric(joined$HRD_Score)

# Outcomes of interest (to predict)
# HRD_Score
# Extract relevant columns
joined.relevant <- data.frame(joined[,c(1, 4:610)], 
                              joined$HRD_Score)
joined.relevant <- joined.relevant %>% dplyr::rename(HRD_Score = joined.HRD_Score)
dim(joined.relevant)

# Throw out rows with high missigness, otherwise mean imputation
na.rows <- joined.relevant[rowSums(is.na(joined.relevant)) > 0,] # two rows comprising all NA values
joined.cleaned <- joined.relevant[-as.numeric(rownames(na.rows)),]
dim(joined.cleaned)
sum(is.na(joined.cleaned)) # verify 0 missingness

# Plot HRD score distribution with GMM (gaussian mixture model)
gmm.fit <- Mclust(joined.cleaned$HRD_Score, G=4, model="V")
summary(gmm.fit, parameters=TRUE)
# plot(gmm.fit, what="density", main="", xlab="HRD Score")
joined.cleaned$GMM_Cluster <- gmm.fit$classification
joined.cleaned$GMM_Cluster <- mapvalues(joined.cleaned$GMM_Cluster,
                                        from=c(1,2,3,4),
                                        to=c('r', 'r', 'r', 'b'))

# Manually find intersection between two density plots
# Plot HRD score distribution with chosen threshold
HRD_Score.threshold <- 44.5
p <- ggplot(joined.cleaned[joined.cleaned$HRD_Score >= 0, ], aes(x=HRD_Score)) + geom_histogram(aes(y=..density..), color='black', fill='darkgrey', bins=200)
p + stat_density(aes(x=HRD_Score), alpha=0.5) + geom_density(aes(x=HRD_Score, color=GMM_Cluster, fill=GMM_Cluster), alpha = 0.5) + theme(legend.position = "none") + geom_vline(aes(xintercept=HRD_Score.threshold), color="red", linetype="dashed", size=1)
# p + geom_density(alpha=0.5, fill='grey') + geom_vline(aes(xintercept=HRD_Score.threshold), color="red", linetype="dashed", size=1)
sum(joined.cleaned$HRD_Score > HRD_Score.threshold)/length(joined.cleaned$HRD_Score) # Percent above threshold cutoff

# After plotting, remove GMM_Cluster column
joined.cleaned = joined.cleaned[,-610]

# Convert continuous (regression) label --> binary (classification) label
# Above / below Z-score = 0
joined.cleaned$HRD_Score_Binarized <- ifelse(joined.cleaned$HRD_Score > HRD_Score.threshold, 1, 0)

# Extract Tissue Source Site Info (TSS)
joined.cleaned$TSS <- substr(joined.cleaned$bcr_patient_barcode, 6, 7)
site.breakdown <- as.data.frame(table(joined.cleaned$TSS))
colnames(site.breakdown) <- c('TSS', 'Sample Count')

# Hold-out TSS's
tss.holdout <- c('BH', 'A2')
use <- dplyr::filter(joined.cleaned, !joined.cleaned$TSS %in% tss.holdout)
holdout <- dplyr::filter(joined.cleaned, joined.cleaned$TSS %in% tss.holdout)
dim(use)[1]/dim(joined.cleaned)[1]
dim(holdout)[1]/dim(joined.cleaned)[1]

sum(joined.cleaned$HRD_Score_Binarized)/dim(joined.cleaned)[1]
sum(use$HRD_Score_Binarized)/dim(use)[1]
sum(holdout$HRD_Score_Binarized)/dim(holdout)[1]

# Split into hold-out and to-be-used datasets
# Split defined PER OUTCOME to ensure class-balance
# Save split dataframes
# Split defined PER OUTCOME to ensure class-balance
# Save split dataframes
hrd.use <- use[,-c(611)]
hrd.hold.out <- holdout[,-c(611)]
saveRDS(hrd.use, file="data/genetic_data/HRD/datasets/brca.HRD_Score.joined.USE.rds")
write.csv(hrd.use, file="data/genetic_data/HRD/datasets/brca.HRD_Score.joined.USE.csv")
saveRDS(hrd.hold.out, file="data/genetic_data/HRD/datasets/brca.HRD_Score.joined.HOLDOUT.rds")
write.csv(hrd.hold.out, file="data/genetic_data/HRD/datasets/brca.HRD_Score.joined.HOLDOUT.csv")

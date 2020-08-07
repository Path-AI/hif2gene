library(devtools)
library(readxl)  
library(cluster)
library(stringr)
library(dplyr)

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

# Convert unit of analysis from sample-level to patient-level
# Choose one sample / slide to represent each patient based on largest tumor area
brca.barcodes <- brca.hifs$bcr_patient_barcode 
length(brca.barcodes[duplicated(brca.barcodes)]) # Duplicated barcodes: same patient, different samples
brca.hifs <- brca.hifs %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(brca.hifs)
brca.hifs <- brca.hifs %>% group_by(bcr_patient_barcode) %>% sample_n(1) # Second filter, if tumor areas are equal, choose randomly
dim(brca.hifs)

skcm.barcodes <- skcm.hifs$bcr_patient_barcode 
length(skcm.barcodes[duplicated(skcm.barcodes)]) # Duplicated barcodes: same patient, different samples
skcm.hifs <- skcm.hifs %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`))
dim(skcm.hifs)
skcm.hifs <- skcm.hifs %>% group_by(bcr_patient_barcode) %>% sample_n(1) # Second filter, if tumor areas are equal, choose randomly
dim(skcm.hifs)

stad.barcodes <- stad.hifs$bcr_patient_barcode 
length(stad.barcodes[duplicated(stad.barcodes)]) # Duplicated barcodes: same patient, different samples
stad.hifs <- stad.hifs %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(stad.hifs)
stad.hifs <- stad.hifs %>% group_by(bcr_patient_barcode) %>% sample_n(1) # Second filter, if tumor areas are equal, choose randomly
dim(stad.hifs)

lusc.barcodes <- lusc.hifs$bcr_patient_barcode 
length(lusc.barcodes[duplicated(lusc.barcodes)]) # Duplicated barcodes: same patient, different samples
lusc.hifs <- lusc.hifs %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_FFPE`))
dim(lusc.hifs)
lusc.hifs <- lusc.hifs %>% group_by(bcr_patient_barcode) %>% sample_n(1) # Second filter, if tumor areas are equal, choose randomly
dim(lusc.hifs)

luad.barcodes <- luad.hifs$bcr_patient_barcode 
length(luad.barcodes[duplicated(luad.barcodes)]) # Duplicated barcodes: same patient, different samples
luad.hifs <- luad.hifs %>% group_by(bcr_patient_barcode) %>% filter(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`==max(`AREA (MM2) OF [TUMOR] IN [TISSUE]_HE`))
dim(luad.hifs)
luad.hifs <- luad.hifs %>% group_by(bcr_patient_barcode) %>% sample_n(1) # Second filter, if tumor areas are equal, choose randomly
dim(luad.hifs)

# Note that columns are ordered different, but luckily, rbind will match based on column names
sum(colnames(brca.hifs) == colnames(stad.hifs))
sum(colnames(brca.hifs) == colnames(skcm.hifs))

# Convert PPFE --> HE suffix
colnames(lusc.hifs) <- lapply(colnames(lusc.hifs), str_replace, "_FFPE", "_HE") 
colnames(lusc.hifs)[2] <- 'H & E_ID'
colnames(skcm.hifs) <- lapply(colnames(skcm.hifs), str_replace, "_FFPE", "_HE") 
colnames(skcm.hifs)[2] <- 'H & E_ID'

# Append all cancer subtypes into pan-cancer dataframe
pancancer.hifs <- rbind(brca.hifs, stad.hifs, luad.hifs, lusc.hifs, skcm.hifs)

# Extract all HIF columns + keep cancer subtype column
pancancer.hif.cols <- as.data.frame(pancancer.hifs[,c(3:609, 617)])

# Throw out rows with missigness
na.rows <- pancancer.hif.cols[rowSums(is.na(pancancer.hif.cols)) > 0,]
pancancer.hif.cols <- pancancer.hif.cols[-as.numeric(rownames(na.rows)),]
dim(pancancer.hif.cols)
sum(is.na(pancancer.hif.cols)) # verify 0 missingness
# write.csv(pancancer.hif.cols, file="data/hif_clusters/pancancer.hifs.csv")

#####################
# Generate Clusters #
#####################
h.cutoff <- 0.95
# num.clusters.desired <- 20 # k

# PANCANCER
# Compute column correlation coefficients in form of matrix
pancancer.hif.corr.mat <- cor(pancancer.hif.cols[,-608], method = c("spearman"))
dim(pancancer.hif.corr.mat)
# Compute correlation distance 
pancancer.hif.corr.dist <- as.dist(1-abs(pancancer.hif.corr.mat))

# Agglomerative clustering
agg.clust <- agnes(pancancer.hif.corr.dist, method = "complete") %>% as.hclust 
plot(agg.clust, which.plots=2, cex = 0.1, hang = -1)
# Define sparse group lasso (SGL) group membership based on h threshold / cutoff
groups <- cutree(agg.clust, h=h.cutoff)
# groups <- cutree(agg.clust, k=num.clusters.desired)
num.groups <- groups %>% unique %>% length
out <- sapply(1:num.groups, function(i) names(groups[groups == i]) %>% print(quote = F))
table(groups)
saveRDS(groups, file="data/hif_clusters/pancancer.hif.clusters.rds") # output R datafile
write.csv(groups, file="data/hif_clusters/pancancer.hif.clusters.csv") # output CSV

# Draw rectangles onto agg.clust plot to denote groupings
rect.hclust(agg.clust, k = num.groups, border = 2:(num.groups+1)) 

# Plot number clusters vs. max intercluster correlation
plot(80:99/100, sapply(80:99, function(i) cutree(agg.clust %>% as.hclust, h=i/100) 
                       %>% table %>% length), xlab = "h (1 - Max Intercluster Correlation)", ylab = "Cluster Count", xlim = c(0.8, 1.00), type = 'o')

# CANCER SUBTYPES
cancer.subtype <- 'BRCA'
cancer.subtype.hif.cols <- pancancer.hif.cols[pancancer.hif.cols$type==cancer.subtype,-608]
dim(cancer.subtype.hif.cols)

# Compute column correlation coefficients in form of matrix
cancer.subtype.hif.corr.mat <- cor(cancer.subtype.hif.cols, method = c("spearman"))
dim(cancer.subtype.hif.corr.mat)
# Compute correlation distance 
cancer.subtype.hif.corr.dist <- as.dist(1-abs(cancer.subtype.hif.corr.mat))

# Agglomerative clustering
agg.clust <- agnes(cancer.subtype.hif.corr.dist, method = "complete") %>% as.hclust 
plot(agg.clust, which.plots=2, cex = 0.1, hang = -1)
# Define sparse group lasso (SGL) group membership based on h threshold / cutoff
groups <- cutree(agg.clust, h=h.cutoff)
# groups <- cutree(agg.clust, k=num.clusters.desired)
num.groups <- groups %>% unique %>% length
out <- sapply(1:num.groups, function(i) names(groups[groups == i]) %>% print(quote = F))
table(groups)
saveRDS(groups, file=sprintf("data/hif_clusters/%s.hif.clusters.rds", cancer.subtype)) # output R datafile
write.csv(groups, file=sprintf("data/hif_clusters/%s.hif.clusters.csv", cancer.subtype)) # output CSV

# Draw rectangles onto agg.clust plot to denote groupings
rect.hclust(agg.clust, k = num.groups, border = 2:(num.groups+1)) 

# Plot number clusters vs. max intercluster correlation
plot(80:99/100, sapply(80:99, function(i) cutree(agg.clust %>% as.hclust, h=i/100) 
                       %>% table %>% length), xlab = "h (1 - Max Intercluster Correlation)", ylab = "Cluster Count", xlim = c(0.8, 1.00), type = 'o')


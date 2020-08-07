library(devtools)
library(readxl)  
library(cluster)
library(stringr)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

# Define outcome of interest
outcome.of.interest <- 'HRD_Score'

#############
# TRAIN SET #
#############

# Load in joined HIF + DDR dataframes, allocated for training (e.g. not hold-out)
# Already cleaned / converted to patient-level / missingness-removed
brca.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/brca.%s.joined.USE.rds", outcome.of.interest)))
stad.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/stad.%s.joined.USE.rds", outcome.of.interest)))
lusc.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/lusc.%s.joined.USE.rds", outcome.of.interest)))
luad.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/luad.%s.joined.USE.rds", outcome.of.interest)))

dim(brca.df)
dim(stad.df)
dim(lusc.df)
dim(luad.df)
# dim(skcm.df)
# total.samples = dim(brca.df)[1] + dim(stad.df)[1] + dim(lusc.df)[1] + dim(luad.df)[1] + dim(skcm.df)[1]
total.samples = dim(brca.df)[1] + dim(stad.df)[1] + dim(lusc.df)[1] + dim(luad.df)[1]
print(total.samples)

# Note that columns are ordered different, but luckily, rbind will match based on column names
sum(colnames(brca.df) == colnames(lusc.df))

# Convert PPFE --> HE suffix
colnames(lusc.df) <- lapply(colnames(lusc.df), str_replace, "_FFPE", "_HE") 
# colnames(skcm.df) <- lapply(colnames(skcm.df), str_replace, "_FFPE", "_HE") 

# Append all cancer subtypes into pan-cancer dataframe
# pancancer.df <- rbind(brca.df, stad.df, luad.df, lusc.df, skcm.df)
pancancer.df <- rbind(brca.df, stad.df, luad.df, lusc.df)
saveRDS(pancancer.df, file=sprintf("data/genetic_data/HRD/datasets/pancancer.%s.joined.USE.rds", outcome.of.interest))
write.csv(pancancer.df, file=sprintf("data/genetic_data/HRD/datasets/pancancer.%s.joined.USE.csv", outcome.of.interest))

# Summary Stats
sum(pancancer.df$HRD_Score_Binarized) / dim(pancancer.df)[1]

################
# HOLD-OUT SET #
################
# Load in joined HIF + DDR dataframes, allocated for training (e.g. not hold-out)
# Already cleaned / converted to patient-level / missingness-removed
brca.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/brca.%s.joined.HOLDOUT.rds", outcome.of.interest)))
stad.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/stad.%s.joined.HOLDOUT.rds", outcome.of.interest)))
lusc.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/lusc.%s.joined.HOLDOUT.rds", outcome.of.interest)))
luad.df <- as.data.frame(readRDS(file=sprintf("data/genetic_data/HRD/datasets/luad.%s.joined.HOLDOUT.rds", outcome.of.interest)))

dim(brca.df)
dim(stad.df)
dim(lusc.df)
dim(luad.df)
# dim(skcm.df)
# total.samples = dim(brca.df)[1] + dim(stad.df)[1] + dim(lusc.df)[1] + dim(luad.df)[1] + dim(skcm.df)[1]
total.samples = dim(brca.df)[1] + dim(stad.df)[1] + dim(lusc.df)[1] + dim(luad.df)[1]
print(total.samples)

# Note that columns are ordered different, but luckily, rbind will match based on column names
sum(colnames(brca.df) == colnames(lusc.df))

# Convert PPFE --> HE suffix
colnames(lusc.df) <- lapply(colnames(lusc.df), str_replace, "_FFPE", "_HE") 
# colnames(skcm.df) <- lapply(colnames(skcm.df), str_replace, "_FFPE", "_HE") 

# Append all cancer subtypes into pan-cancer dataframe
# pancancer.df <- rbind(brca.df, stad.df, luad.df, lusc.df, skcm.df)
pancancer.df <- rbind(brca.df, stad.df, luad.df, lusc.df)
saveRDS(pancancer.df, file=sprintf("data/genetic_data/HRD/datasets/pancancer.%s.joined.HOLDOUT.rds", outcome.of.interest))
write.csv(pancancer.df, file=sprintf("data/genetic_data/HRD/datasets/pancancer.%s.joined.HOLDOUT.csv", outcome.of.interest))

# Summary Stats
sum(pancancer.df$HRD_Score_Binarized) / dim(pancancer.df)[1]


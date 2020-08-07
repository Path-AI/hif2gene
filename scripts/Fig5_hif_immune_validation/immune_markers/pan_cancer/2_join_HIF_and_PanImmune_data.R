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

# Load in HIF dataframes (filtered for samples / patients with PanImmune data)
brca.pan.immune.joined <- readRDS(file="data/genetic_data/PanImmune/processed/brca.pan.immune.joined.rds")
stad.pan.immune.joined <- readRDS(file="data/genetic_data/PanImmune/processed/stad.pan.immune.joined.rds")
luad.pan.immune.joined <- readRDS(file="data/genetic_data/PanImmune/processed/luad.pan.immune.joined.rds")
lusc.pan.immune.joined <- readRDS(file="data/genetic_data/PanImmune/processed/lusc.pan.immune.joined.rds")
skcm.pan.immune.joined <- readRDS(file="data/genetic_data/PanImmune/processed/skcm.pan.immune.joined.rds")

# Note that columns are ordered different, but luckily, rbind will match based on column names
sum(colnames(brca.pan.immune.joined) == colnames(stad.pan.immune.joined))
sum(colnames(brca.pan.immune.joined) == colnames(skcm.pan.immune.joined))

# Convert PPFE --> HE suffix
colnames(lusc.pan.immune.joined) <- lapply(colnames(lusc.pan.immune.joined), str_replace, "_FFPE", "_HE") 
colnames(lusc.pan.immune.joined)[3] <- 'H & E_ID'
colnames(skcm.pan.immune.joined) <- lapply(colnames(skcm.pan.immune.joined), str_replace, "_FFPE", "_HE") 
colnames(skcm.pan.immune.joined)[3] <- 'H & E_ID'

# Append all cancer subtypes into pan-cancer dataframe
pancancer.pan.immune.joined <- rbind(brca.pan.immune.joined, stad.pan.immune.joined,
                                     luad.pan.immune.joined, lusc.pan.immune.joined, 
                                     skcm.pan.immune.joined)
write.csv(pancancer.pan.immune.joined, file="data/genetic_data/PanImmune/processed/pancancer.pan.immune.joined.csv")


library(devtools)
library(readxl)  
library(caret)
library(dplyr)

# Close any active connections and already loaded libraries 
set.seed(12261995)
closeAllConnections()
rm(list = ls())

# Set working directory
setwd("~/Desktop/hif2gene/")

##########
# DEFINE #
##########
outcome.of.interest <- 'TIGIT'

# Load in ensemble model coefficients = mean(model1, model2, model3) beta coefficients
df <- read.table(sprintf("data/model_outputs/tigit/coefficients/skcm/%s_model_ensemble.csv", outcome.of.interest), sep=",", header=TRUE)

# Convert betas to magnitudes / absolute values
df$AbsBeta <- abs(df$Beta)

# Filter out HIFs with beta = 0
df.filtered = df[df$AbsBeta != 0,]

# Group by HIF cluster, compute summary stats on ensemble absolute beta values
df.agg <- df.filtered %>% dplyr::group_by(Cluster) %>% dplyr::summarise(
  IQR_25=quantile(AbsBeta, p=0.25),
  Median=median(AbsBeta, na.rm = TRUE),
  IQR_75=quantile(AbsBeta, p=0.75),
  Mean=mean(AbsBeta),
  Maximum=max(AbsBeta),
  N=n()
)

df.max <- df.filtered %>% group_by(Cluster) %>% top_n(1, AbsBeta)
df.filtered$Cluster = factor(df.filtered$Cluster) # convert to factor

# Box plot (median, IQR)
# Find top N most "predictive" HIF clusters
# Sort by max
N <- 5
df.agg <- df.agg %>%dplyr::arrange(desc(Maximum))
top.n.hif.clusters <- df.agg$Cluster[1:N]
p <- ggplot(df.filtered[df.filtered$Cluster %in% top.n.hif.clusters, ], aes(x=reorder(Cluster, AbsBeta, FUN=max), y=AbsBeta))
p <- p + geom_jitter(color="black", fill="darkgrey", size=1.2, alpha=0.6) + geom_boxplot(outlier.shape=NA, width=0.5, color="black", fill="#E16C6C")
p <- p + theme(axis.text=element_text(size=16), axis.title.x = element_blank(), axis.title.y = element_blank(), 
               axis.title=element_text(size=14)) + coord_flip()
p


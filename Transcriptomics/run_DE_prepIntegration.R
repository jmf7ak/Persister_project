##### load packages #####
library(tidyverse)

##### load data #####
# set working directory
setwd('..')
dir <- getwd()

# load files
T5T0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T5T0_results.csv"))
U5U0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U5U0_results.csv"))
T24T0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T24T0_results.csv"))
U24U0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U24U0_results.csv"))

# find model directory
setwd('..')
model_dir <- getwd()
model_dir <- paste0(model_dir, '/', 'Model')

metabGenes <- read.csv(file.path(model_dir, 'data', 'PA14recon1_v25_genes.csv'))

# filter all DE data with for metab genes
T5T0_all2metab <- T5T0[which(T5T0$X %in% metabGenes$PA14_metabolic_model_genes),]
T24T0_all2metab <- T24T0[which(T24T0$X %in% metabGenes$PA14_metabolic_model_genes),]
U5U0_all2metab <- U5U0[which(U5U0$X %in% metabGenes$PA14_metabolic_model_genes),]
U24U0_all2metab <- U24U0[which(U24U0$X %in% metabGenes$PA14_metabolic_model_genes),]

# re-order according to padj and log2FoldChange
T5T0_all2metab <- T5T0_all2metab[order(T5T0_all2metab$padj),]
T24T0_all2metab <- T24T0_all2metab[order(T24T0_all2metab$padj),]
U5U0_all2metab <- U5U0_all2metab[order(U5U0_all2metab$padj),]
U24U0_all2metab <- U24U0_all2metab[order(U24U0_all2metab$padj),]

# count how many genes have a padj < 0.5
sum(T5T0_all2metab$padj < 0.5, na.rm=TRUE)
sum(T24T0_all2metab$padj < 0.5, na.rm=TRUE)
sum(U5U0_all2metab$padj < 0.5, na.rm=TRUE)
sum(U24U0_all2metab$padj < 0.5, na.rm=TRUE)

# count how many genes have a padj < 0.01
sum(T5T0_all2metab$padj < 0.01, na.rm=TRUE)
sum(T24T0_all2metab$padj < 0.01, na.rm=TRUE)
sum(U5U0_all2metab$padj < 0.01, na.rm=TRUE)
sum(U24U0_all2metab$padj < 0.01, na.rm=TRUE)

# count how many genes have a abs(log2FoldChange) > 2 and padj < 0.01
sum(abs(T5T0_all2metab$log2FoldChange) > 2 & T5T0_all2metab$padj < 0.01, na.rm=TRUE)
sum(abs(T24T0_all2metab$log2FoldChange) > 2 & T24T0_all2metab$padj < 0.01, na.rm=TRUE)
sum(abs(U5U0_all2metab$log2FoldChange) > 2 & U5U0_all2metab$padj < 0.01, na.rm=TRUE)
sum(abs(U24U0_all2metab$log2FoldChange) > 2 & U24U0_all2metab$padj < 0.01, na.rm=TRUE)

# count how many genes have padj == NA
sum(is.na(T5T0_all2metab$padj))
sum(is.na(T24T0_all2metab$padj))
sum(is.na(U5U0_all2metab$padj))
sum(is.na(U24U0_all2metab$padj))

# create dataframes for integration in PA14 models with removed NAs
T5T0_all2metab <- T5T0_all2metab[order(T5T0_all2metab$X),]
T24T0_all2metab <- T24T0_all2metab[order(T24T0_all2metab$X),]
U5U0_all2metab <- U5U0_all2metab[order(U5U0_all2metab$X),]
U24U0_all2metab <- U24U0_all2metab[order(U24U0_all2metab$X),]

all2metab_integration <- data.frame(genes = T5T0_all2metab[complete.cases(T5T0_all2metab),]$X,
                                    T5T0_all2metab.log2FoldChange = T5T0_all2metab[complete.cases(T5T0_all2metab),]$log2FoldChange, 
                                    T5T0_all2metab.padj = T5T0_all2metab[complete.cases(T5T0_all2metab),]$padj,
                                    T24T0_all2metab.log2FoldChange = T24T0_all2metab[complete.cases(T24T0_all2metab),]$log2FoldChange, 
                                    T24T0_all2metab.padj = T24T0_all2metab[complete.cases(T24T0_all2metab),]$padj,
                                    U5U0_all2metab.log2FoldChange = U5U0_all2metab[complete.cases(U5U0_all2metab),]$log2FoldChange, 
                                    U5U0_all2metab.padj = U5U0_all2metab[complete.cases(U5U0_all2metab),]$padj,
                                    U24U0_all2metab.log2FoldChange = U24U0_all2metab[complete.cases(U24U0_all2metab),]$log2FoldChange, 
                                    U24U0_all2metab.padj = U24U0_all2metab[complete.cases(U24U0_all2metab),]$padj)

write.csv(as.data.frame(all2metab_integration),file=file.path(model_dir, 'data', "DESeq_dds_groupbatch_all2metab_integration_180221.csv"))

# create dataframes for integration in PA14 models with NAs
T5T0_all2metab <- T5T0_all2metab[order(T5T0_all2metab$X),]
T24T0_all2metab <- T24T0_all2metab[order(T24T0_all2metab$X),]
U5U0_all2metab <- U5U0_all2metab[order(U5U0_all2metab$X),]
U24U0_all2metab <- U24U0_all2metab[order(U24U0_all2metab$X),]

all2metab_integration_NA <- data.frame(genes = T5T0_all2metab$X,
                                       T5T0_all2metab.log2FoldChange = T5T0_all2metab$log2FoldChange, 
                                       T5T0_all2metab.padj = T5T0_all2metab$padj,
                                       T24T0_all2metab.log2FoldChange = T24T0_all2metab$log2FoldChange, 
                                       T24T0_all2metab.padj = T24T0_all2metab$padj,
                                       U5U0_all2metab.log2FoldChange = U5U0_all2metab$log2FoldChange, 
                                       U5U0_all2metab.padj = U5U0_all2metab$padj,
                                       U24U0_all2metab.log2FoldChange = U24U0_all2metab$log2FoldChange, 
                                       U24U0_all2metab.padj = U24U0_all2metab$padj)

write.csv(as.data.frame(all2metab_integration_NA),file=file.path(model_dir, 'data', "DESeq_dds_groupbatch_all2metab_integration_NA_180221.csv"))

# basic stats
sum(T5T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(T24T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U5U0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U24U0_all2metab$padj < 0.05, na.rm=TRUE)


sum(T5T0_all2metab$log2FoldChange > 0 & T5T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(T24T0_all2metab$log2FoldChange > 0 & T24T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U5U0_all2metab$log2FoldChange > 0 & U5U0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U24U0_all2metab$log2FoldChange > 0 & U24U0_all2metab$padj < 0.05, na.rm=TRUE)

sum(T5T0_all2metab$log2FoldChange < 0 & T5T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(T24T0_all2metab$log2FoldChange < 0 & T24T0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U5U0_all2metab$log2FoldChange < 0 & U5U0_all2metab$padj < 0.05, na.rm=TRUE)
sum(U24U0_all2metab$log2FoldChange < 0 & U24U0_all2metab$padj < 0.05, na.rm=TRUE)
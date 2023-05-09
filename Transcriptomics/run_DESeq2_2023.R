##### load packages #####
library(DESeq2)
library(vsn)
library(pheatmap)
library(grid)
library(RColorBrewer)

##### load data #####
# set working directory
setwd('C:/Users/jmfic/OneDrive - University of Virginia/21_SUMMA/persister/Transcriptomics')
dir <- getwd()

# obtain names of all files that include ".counts" in the name
files <- list.files(path = file.path(dir, 'analysis', 'counts2'), pattern="*.counts")

# provide the names of all the files *note that sampleName must be in the same order as files
sampleName <- c("treated_0A_81017d",
                "treated_0B_81017d",
                "treated_0C_81017d",
                "treated_0F_81017d",
                "untreated_0A_81017d",
                "untreated_0B_81017d",
                "untreated_0C_81017d",
                "untreated_0D_81017d", 
                
                
                "treated_24A_81017d",
                "treated_24B_81017d",
                "treated_24C_81017d",
                "treated_24F_81017d",
                "untreated_24A_81017d",
                "untreated_24B_81017d",
                "untreated_24C_81017d",
                "untreated_24D_81017d",
                
                "treated_5A_81017d",
                "treated_5B_81017d",
                "treated_5C_81017d",
                "treated_5F_81017d",
                "untreated_5A_81017d",
                "untreated_5B_81017d",
                "untreated_5C_81017d",
                "untreated_5D_81017d",
                
                "treated_0A_83017d",
                "treated_0B_83017d",
                "treated_0D_83017d",
                "treated_0F_83017d",
                "untreated_0A_83017d",
                "untreated_0B_83017d",
                #"untreated_0C_83017d", ###corrupt
               # "untreated_0D_83017d", ###corrupt
                
               #"treated_24A_83017d",
               "treated_24B_83017d",
               "treated_24D_83017d",
               "treated_24F_83017d",
               "untreated_24A_83017d",
               "untreated_24B_83017d",
               "untreated_24C_83017d",
               "untreated_24D_83017d",
               
                "treated_5A_83017d",
                "treated_5B_83017d",
                "treated_5D_83017d",
                "treated_5F_83017d",
                "untreated_5A_83017d",
                "untreated_5B_83017d",
                "untreated_5C_83017d",
                "untreated_5D_83017d")
                
               

# load in files as a table and assign them their sampleName
for (i in 1:length(files)) assign(sampleName[i], read.table(file.path(dir, 'analysis', 'counts2', files[i])))

##### dds set-up #####
# create a gene_id dataframe
gene_id <- data.frame(gene_id = treated_0A_81017d$V1)
# remove the first row
gene_id <- gene_id[-1,]
# create a counts dataframe
counts <- data.frame(treated_0A_81017d = treated_0A_81017d$V7,
                     treated_0B_81017d = treated_0B_81017d$V7,
                     treated_0C_81017d = treated_0C_81017d$V7,
                     treated_0F_81017d = treated_0F_81017d$V7,
                     untreated_0A_81017d = untreated_0A_81017d$V7,
                     untreated_0B_81017d = untreated_0B_81017d$V7,
                     untreated_0C_81017d = untreated_0C_81017d$V7,
                     untreated_0D_81017d = untreated_0D_81017d$V7,
                     
                     treated_5A_81017d = treated_5A_81017d$V7,
                     treated_5B_81017d = treated_5B_81017d$V7,
                     treated_5C_81017d = treated_5C_81017d$V7,
                     treated_5F_81017d = treated_5F_81017d$V7,
                     untreated_5A_81017d = untreated_5A_81017d$V7,
                     untreated_5B_81017d = untreated_5B_81017d$V7,
                     untreated_5C_81017d = untreated_5C_81017d$V7,
                     untreated_5D_81017d = untreated_5D_81017d$V7,
                     
                     treated_24A_81017d = treated_24A_81017d$V7,
                     treated_24B_81017d = treated_24B_81017d$V7,
                     treated_24C_81017d = treated_24C_81017d$V7,
                     treated_24F_81017d = treated_24F_81017d$V7,
                     untreated_24A_81017d = untreated_24A_81017d$V7,
                     untreated_24B_81017d = untreated_24B_81017d$V7,
                     untreated_24C_81017d = untreated_24C_81017d$V7,
                     untreated_24D_81017d = untreated_24D_81017d$V7,
                     
                     treated_0A_83017d = treated_0A_83017d$V7,
                     treated_0B_83017d = treated_0B_83017d$V7,
                     treated_0D_83017d = treated_0D_83017d$V7,
                     treated_0F_83017d = treated_0F_83017d$V7,
                     untreated_0A_83017d = untreated_0A_83017d$V7,
                     untreated_0B_83017d = untreated_0B_83017d$V7,
                     #untreated_0C_83017d = untreated_0C_83017d$V7,
                     #untreated_0D_83017d = untreated_0D_83017d$V7,
                     
                     treated_5A_83017d = treated_5A_83017d$V7,
                     treated_5B_83017d = treated_5B_83017d$V7,
                     treated_5D_83017d = treated_5D_83017d$V7,
                     treated_5F_83017d = treated_5F_83017d$V7,
                     untreated_5A_83017d = untreated_5A_83017d$V7,
                     untreated_5B_83017d = untreated_5B_83017d$V7,
                     untreated_5C_83017d = untreated_5C_83017d$V7,
                     untreated_5D_83017d = untreated_5D_83017d$V7,
                     
                     #treated_24A_83017d = treated_24A_83017d$V7,
                     treated_24B_83017d = treated_24B_83017d$V7,
                     treated_24D_83017d = treated_24D_83017d$V7,
                     treated_24F_83017d = treated_24F_83017d$V7,
                     untreated_24A_83017d = untreated_24A_83017d$V7,
                     untreated_24B_83017d = untreated_24B_83017d$V7,
                     untreated_24C_83017d = untreated_24C_83017d$V7,
                     untreated_24D_83017d = untreated_24D_83017d$V7)
# remove the first row
counts <- counts[-1,]
# add gene_ids to rownames
rownames(counts) <- gene_id
# convert to a matrix so that the innards are numeric
counts <- data.matrix(counts)
# create a new dataframe, coldata, that contains condition information
coldata <- data.frame(condition = c("treated","treated","treated","treated",
                                    "untreated","untreated","untreated","untreated",
                                    "treated","treated","treated","treated",
                                    "untreated","untreated","untreated","untreated",
                                    "treated","treated","treated","treated",
                                    "untreated","untreated","untreated","untreated",
                                    "treated","treated","treated","treated",
                                    "untreated","untreated",
                                    "treated","treated","treated","treated",
                                    "untreated","untreated","untreated","untreated",
                                    "treated","treated","treated",
                                    "untreated","untreated","untreated","untreated"),
                      time = c("0","0","0","0",
                               "0","0","0","0",
                               "5","5","5","5",
                               "5","5","5","5",
                               "24","24","24","24",
                               "24","24","24","24",
                               "0","0","0","0",
                               "0","0",
                               "5","5","5","5",
                               "5","5","5","5",
                               "24","24","24",
                               "24","24","24","24"),
                      batch = c("081017","081017","081017","081017",
                                "081017","081017","081017","081017",
                                "081017","081017","081017","081017",
                                "081017","081017","081017","081017",
                                "081017","081017","081017","081017",
                                "081017","081017","081017","081017",
                                "083017","083017","083017","083017",
                                "083017","083017",
                                "083017","083017","083017","083017",
                                "083017","083017","083017","083017",
                                "083017","083017","083017",
                                "083017","083017","083017","083017"))
# rename the rows of coldata such that they match the sampleName
rownames(coldata) = c("treated_0A_81017d",
                      "treated_0B_81017d",
                      "treated_0C_81017d",
                      "treated_0F_81017d",
                      "untreated_0A_81017d",
                      "untreated_0B_81017d",
                      "untreated_0C_81017d",
                      "untreated_0D_81017d",
                      
                      "treated_5A_81017d",
                      "treated_5B_81017d",
                      "treated_5C_81017d",
                      "treated_5F_81017d",
                      "untreated_5A_81017d",
                      "untreated_5B_81017d",
                      "untreated_5C_81017d",
                      "untreated_5D_81017d",
                      
                      "treated_24A_81017d",
                      "treated_24B_81017d",
                      "treated_24C_81017d",
                      "treated_24F_81017d",
                      "untreated_24A_81017d",
                      "untreated_24B_81017d",
                      "untreated_24C_81017d",
                      "untreated_24D_81017d",
                      
                      "treated_0A_83017d",
                      "treated_0B_83017d",
                      "treated_0D_83017d",
                      "treated_0F_83017d",
                      "untreated_0A_83017d",
                      "untreated_0B_83017d",
                    #  "untreated_0C_83017d",
                      #"untreated_0D_83017d",
                      
                      "treated_5A_83017d",
                      "treated_5B_83017d",
                      "treated_5D_83017d",
                      "treated_5F_83017d",
                      "untreated_5A_83017d",
                      "untreated_5B_83017d",
                      "untreated_5C_83017d",
                      "untreated_5D_83017d",
                      
                      #"treated_24A_83017d",
                      "treated_24B_83017d",
                      "treated_24D_83017d",
                      "treated_24F_83017d",
                      "untreated_24A_83017d",
                      "untreated_24B_83017d",
                      "untreated_24C_83017d",
                      "untreated_24D_83017d")
# check to see if rownames of coldata match colnames of counts
all(rownames(coldata) %in% colnames(counts))
# filter the counts dataframe so that it only contains the columns that match the rownames of coldata *note, this should not be necessary ...
# based on the above commands
counts <- counts[, rownames(coldata)]
# check to see if rownames of coldata are in the same order of colnames of counts
all(rownames(coldata) == colnames(counts))
# for multi-factorial design, join the two factors of interest into a "group" variable
coldata$group <- factor(paste0(coldata$time,coldata$condition))

##### design ~ batch + condition #####
# construct a DESeqDataSet, dds, using counts and coldata, and set the design to compare groups controlling for batch effects
dds_batchgroup <- DESeqDataSetFromMatrix(countData = counts,
                                          colData = coldata,
                                          design = ~ batch + group)
dds_batchgroup

mcols(dds_batchgroup) <- DataFrame(mcols(dds_batchgroup),gene_id)
mcols(dds_batchgroup)

dds_batchgroup <- DESeq(dds_batchgroup)

# obtain comparison for 5treated vs 0treated (5treated/0treated)
T5T0batch <- results(dds_batchgroup, contrast=c("group","5treated","0treated"))
# obtain comparison for 5untreated vs 0untreated (5untreated/0untreated)
U5U0batch <- results(dds_batchgroup, contrast=c("group","5untreated","0untreated"))
# obtain comparison for 24treated vs 0treated (24treated/0treated)
T24T0batch <- results(dds_batchgroup, contrast=c("group","24treated","0treated"))
# obtain comparison for 24untreated vs 0untreated (24untreated/0untreated)
U24U0batch <- results(dds_batchgroup, contrast=c("group","24untreated","0untreated"))
# obtain comparison for 5treated vs 5untreated (5treated/5untreated)
T5U5batch <- results(dds_batchgroup, contrast=c("group","5treated","5untreated"))
# obtain comparison for 24treated vs 24untreated (24treated/24untreated)
T24U24batch <- results(dds_batchgroup, contrast=c("group","24treated","24untreated"))

# get summary information on the differential expression comparison
summary(T5T0batch)
summary(U5U0batch)
summary(T24T0batch)
summary(U24U0batch)
summary(T5U5batch)
summary(T24U24batch)

sum(T5T0batch$padj < 0.05, na.rm=TRUE)
sum(U5U0batch$padj < 0.05, na.rm=TRUE)
sum(T24T0batch$padj < 0.05, na.rm=TRUE)
sum(U24U0batch$padj < 0.05, na.rm=TRUE)
sum(T5U5batch$padj < 0.05, na.rm=TRUE)
sum(T24U24batch$padj < 0.05, na.rm=TRUE)

# MA plot
plotMA(T5T0batch, ylim = c(-4,4), main='T5T0')
plotMA(U5U0batch, ylim = c(-4,4), main='U5U0')
plotMA(T24T0batch, ylim = c(-4,4), main='T24T0')
plotMA(U24U0batch, ylim = c(-4,4), main='U24U0')
plotMA(T5U5batch, ylim = c(-4,4), main='T5U5')
plotMA(T24U24batch, ylim = c(-4,4), main='T24U24')

write.csv(as.data.frame(T5T0batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T5T0_results_joe.csv"))
write.csv(as.data.frame(U5U0batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U5U0_results_joe.csv"))
write.csv(as.data.frame(T24T0batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T24T0_results_joe.csv"))
write.csv(as.data.frame(U24U0batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U24U0_results_joe.csv"))
write.csv(as.data.frame(T5U5batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T5U5_results_joe.csv"))
write.csv(as.data.frame(T24U24batch),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T24U24_results_joe.csv"))

##### transformations ~ data #####

vsd_group <- varianceStabilizingTransformation(dds_batchgroup, blind = FALSE)

notAllZero_group <- (rowSums(counts(dds_batchgroup))>0)
meanSdPlot(assay(vsd_group[notAllZero_group,]),ylim=c(0,5),ylab="vsd - sd")

# heatmap of the count matrix
select_group <- order(rowMeans(counts(dds_batchgroup,normalized=TRUE)),decreasing=FALSE)[1:100]
df_group <- as.data.frame(colData(dds_batchgroup)[,c("group")])
rownames(df_group) <- colnames(dds_batchgroup)
pheatmap(assay(vsd_group)[select_group,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation_col=df_group)

# sample-to-sample distances
sampleDists_vsd <- dist(t(assay(vsd_group)))
sampleDists_vsd <- dist(t(assay(vsd_group)))
sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)
rownames(sampleDistMatrix_vsd) <- rownames(coldata)
colnames(sampleDistMatrix_vsd) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix_vsd,
         clustering_distance_rows=sampleDists_vsd,
         clustering_distance_cols=sampleDists_vsd,
         col=colors)

# principal component plot of the samples
plotPCA(vsd_group, intgroup=c("condition","time","batch"))

##### design ~ time + batch #####

# construct a DESeqDataSet, dds, using counts and coldata, and set the design to include an interaction term between treatment and time
dds_timebatch <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition + time + batch + condition:time)
dds_timebatch$condition <- factor(dds_timebatch$condition, levels = c("untreated","treated"))
dds_timebatch

mcols(dds_timebatch) <- DataFrame(mcols(dds_timebatch),gene_id)
mcols(dds_timebatch)

dds_timebatch<- DESeq(dds_timebatch, test="LRT", reduced = ~ condition + time + batch)

timebatch <- results(dds_timebatch)
summary(timebatch)
sum(timebatch$padj < 0.05, na.rm = TRUE)

# for contrasts, see this forum: https://support.bioconductor.org/p/101002/

TU0 <- results(dds_timebatch, name = "condition_treated_vs_untreated", test="Wald")
sum(TU0$padj < 0.05, na.rm = TRUE)

TU5 <- results(dds_timebatch, contrast = list(c("condition_treated_vs_untreated","conditiontreated.time5")), test="Wald")
sum(TU5$padj < 0.05, na.rm = TRUE)

TU24 <- results(dds_timebatch, contrast = list(c("condition_treated_vs_untreated","conditiontreated.time24")), test="Wald")
sum(TU24$padj < 0.05, na.rm = TRUE)

batch <- results(dds_timebatch, name = "batch_083017_vs_081017", test="Wald")
sum(batch$padj < 0.05, na.rm = TRUE)

write.csv(as.data.frame(TU5),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_timebatch_TU5_results_joe.csv"))
write.csv(as.data.frame(TU24),file=file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_timebatch_TU24_results_joe.csv"))


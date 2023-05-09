##### load packages #####
library(tidyverse)
library(ggdendro)
#library(plotly)
library(gplots)
library(UpSetR)
library(pheatmap)
library(dplyr)

##### create color palette #####
myPalette <- c("#000000","#999999", "#00c590", "#008a65", "#00503A", "#ffcd5d", "#E69F00", "#996900", "#0072B2", "#D55E00","#56B4E9","#CC79A7")

# black, gray, bright green, *green*[4], dark green, light yellow, *yellow*[7], dark yellow, blue, orange, light blue, *pink* [12]
# from here: http://www.colorhexa.com/e79f00
# and here: http://www.colorhexa.com/009e73

##### load data #####

# set working directory
setwd('C:/Users/jmfic/OneDrive - University of Virginia/21_SUMMA/persister/Transcriptomics')
dir <- getwd()

# load files
T5T0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T5T0_results_2023.csv"))
U5U0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U5T0_results_2023.csv"))
T24T0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_T24T0_results_2023.csv"))
U24U0 <- read.csv(file.path(dir, 'analysis', 'DESeq2', "DESeq_dds_groupbatch_U24T0_results_2023.csv"))

T5T0$condition <- c("T5T0")
U5U0$condition <- c("U5U0")
T24T0$condition <- c("T24T0")
U24U0$condition <- c("U24U0")

complete <- rbind(T5T0,U5U0,T24T0,U24U0)
complete$condition <- factor(complete$condition,levels=c("T5T0","U5U0","T24T0","U24U0"))

T5T0_sig <- filter(T5T0,padj < 0.01)
U5U0_sig <- filter(U5U0,padj < 0.01)
T24T0_sig <- filter(T24T0,padj < 0.01)
U24U0_sig <- filter(U24U0,padj < 0.01)

T5T0_sig_up <- filter(T5T0, padj < 0.01 & log2FoldChange > 0)
T5T0_sig_down <- filter(T5T0, padj < 0.01 & log2FoldChange < 0)
U5U0_sig_up <- filter(U5U0, padj < 0.01 & log2FoldChange > 0)
U5U0_sig_down <- filter(U5U0, padj < 0.01 & log2FoldChange < 0)
T24T0_sig_up <- filter(T24T0, padj < 0.01 & log2FoldChange > 0)
T24T0_sig_down <- filter(T24T0, padj < 0.01 & log2FoldChange < 0)
U24U0_sig_up <- filter(U24U0, padj < 0.01 & log2FoldChange > 0)
U24U0_sig_down <- filter(U24U0, padj < 0.01 & log2FoldChange < 0)

write.csv(T5T0_sig_up,file.path(dir, 'analysis', 'DESeq2', 'T5T0_upregulated_2023.csv'))
write.csv(T5T0_sig_down,file.path(dir, 'analysis', 'DESeq2', 'T5T0_downregulated_2023.csv'))
write.csv(U5U0_sig_up,file.path(dir, 'analysis', 'DESeq2', 'U5U0_upregulated_2023.csv'))
write.csv(U5U0_sig_down,file.path(dir, 'analysis', 'DESeq2', 'U5U0_downregulated_2023.csv'))
write.csv(T24T0_sig_up,file.path(dir, 'analysis', 'DESeq2', 'T24T0_upregulated_2023.csv'))
write.csv(T24T0_sig_down,file.path(dir, 'analysis', 'DESeq2', 'T24T0_downregulated_2023.csv'))
write.csv(U24U0_sig_up,file.path(dir, 'analysis', 'DESeq2', 'U24U0_upregulated_2023.csv'))
write.csv(U24U0_sig_down,file.path(dir, 'analysis', 'DESeq2', 'U24U0_downregulated_2023.csv'))

#### visualize overlap ####

upregulated <- list(P5 = T5T0_sig_up$X, P24 = T24T0_sig_up$X, U5 = U5U0_sig_up$X, U24 = U24U0_sig_up$X)
downregulated <- list(P5 = T5T0_sig_down$X, P24 = T24T0_sig_down$X, U5 = U5U0_sig_down$X, U24 = U24U0_sig_down$X)

upset(fromList(upregulated),
      text.scale = c(1.8, 1.9, 1.9, 1.9, 1.9,1.9),
      sets = rev(c('P5','P24','U5','U24')),
      order.by = "freq",
      #empty.intersections = "on",
      point.size = 2, line.size = 0.7,
      show.numbers = "no",
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Upregulated Genes",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("P5","P24"),
                          color = myPalette[4], active = T,
                          query.name = "Unique Persister"
                          ),
                     list(query = intersects,
                          params = list("U5","U24"),
                          color = myPalette[7], active = T,
                          query.name = "Unique Untreated")))

upset(fromList(downregulated),
      text.scale = c(1.8, 1.9, 1.9, 1.9, 1.9,1.9),
      sets = rev(c('P5','P24','U5','U24')),
      order.by = "freq",
      #empty.intersections = "on",
      point.size = 2, line.size = 0.7,
      show.numbers = "no",
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Downregulated Genes",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("P5","P24"),
                          color = myPalette[4], active = T,
                          query.name = "Unique Persister"
                        ),
                     list(query = intersects,
                          params = list("U5","U24"),
                          color = myPalette[7], active = T,
                          query.name = "Unique Untreated")))


# find the list of genes that come up as significantly DE in at least one of the conditions
sig_genes_groupbatch_v0<- unique(select(rbind(T5T0_sig,U5U0_sig,T24T0_sig,U24U0_sig),X))

##### heatmap #####

heatmap_sig_groupbatch_v0 <- sig_genes_groupbatch_v0
heatmap_sig_groupbatch_v0$X <- sort(heatmap_sig_groupbatch_v0$X)
heatmap_sig_groupbatch_v0$T5T0 <- rep(0,nrow(heatmap_sig_groupbatch_v0))
heatmap_sig_groupbatch_v0$U5U0 <- rep(0,nrow(heatmap_sig_groupbatch_v0))
heatmap_sig_groupbatch_v0$T24T0 <- rep(0,nrow(heatmap_sig_groupbatch_v0))
heatmap_sig_groupbatch_v0$U24U0 <- rep(0,nrow(heatmap_sig_groupbatch_v0))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
heatmap_sig_groupbatch_v0[heatmap_sig_groupbatch_v0$X %in% T5T0$X,]$T5T0 = T5T0[T5T0$X %in% heatmap_sig_groupbatch_v0$X,]$log2FoldChange
heatmap_sig_groupbatch_v0[heatmap_sig_groupbatch_v0$X %in% U5U0$X,]$U5U0 = U5U0[U5U0$X %in% heatmap_sig_groupbatch_v0$X,]$log2FoldChange
heatmap_sig_groupbatch_v0[heatmap_sig_groupbatch_v0$X %in% T24T0$X,]$T24T0 = T24T0[T24T0$X %in% heatmap_sig_groupbatch_v0$X,]$log2FoldChange
heatmap_sig_groupbatch_v0[heatmap_sig_groupbatch_v0$X %in% U24U0$X,]$U24U0 = U24U0[U24U0$X %in% heatmap_sig_groupbatch_v0$X,]$log2FoldChange

rownames(heatmap_sig_groupbatch_v0) <- heatmap_sig_groupbatch_v0$X
heatmap_sig_groupbatch_v0 <- as.matrix(select(heatmap_sig_groupbatch_v0,-X))

myBreaks <- c(min(heatmap_sig_groupbatch_v0),seq(-2+1e-5,2-1e-5,length=99),max(heatmap_sig_groupbatch_v0))
hmcols <- c(colorRampPalette(c(myPalette[9],"white",myPalette[10]))(100))

colnames(heatmap_sig_groupbatch_v0) <- c('P5','U5','P24','U24')

pheatmap(heatmap_sig_groupbatch_v0,
         angle_col = "45",
         color = hmcols,
         breaks = myBreaks,
         show_rownames = FALSE,
         treeheight_row = 0,
         cellwidth = 25,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "correlation",
         clustering_method = "average")

##### volcano plot #####
# inspired from here: https://twbattaglia.github.io/2016/12/17/volcano-plot/

complete <- complete %>%
  mutate(color = ifelse(complete$log2FoldChange>2 & complete$padj < 0.01,
                        yes = "num",
                        no = ifelse(complete$log2FoldChange<(-2) & complete$padj < 0.01,
                                    yes = "denom",
                                    no = "none")))

ggplot(complete, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-5,5)) +
  ylim(c(0,50)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

#### enrichment analysis ####

enrich_U5U0_up <- read.csv(file.path(dir, 'analysis', 'enrichment', 'U5U0_up_2023_enrich.tsv.csv'), header = TRUE)
enrich_U5U0_down <- read.csv(file.path(dir, 'analysis', 'enrichment', 'U5U0_down_2023_enrich.tsv.csv'), header = TRUE)
enrich_U24U0_up <- read.csv(file.path(dir, 'analysis', 'enrichment', 'U24U0_up_2023_enrich.tsv.csv'), header = TRUE)
enrich_U24U0_down <- read.csv(file.path(dir, 'analysis', 'enrichment', 'U24U0_down_2023_enrich.tsv.csv'), header = TRUE)
enrich_T5T0_up <- read.csv(file.path(dir, 'analysis', 'enrichment', 'T5T0_up_2023_enrich.tsv.csv'), header = TRUE)
enrich_T5T0_down <- read.csv(file.path(dir, 'analysis', 'enrichment', 'T5T0_down_2023_enrich.tsv.csv'), header = TRUE)
enrich_T24T0_up <- read.csv(file.path(dir, 'analysis', 'enrichment', 'T24T0_up_2023_enrich.tsv.csv'), header = TRUE)
enrich_T24T0_down <- read.csv(file.path(dir, 'analysis', 'enrichment', 'T24T0_down_2023_enrich.tsv.csv'), header = TRUE)

enrich_U5U0_up$perc_set <- (enrich_U5U0_up$genes_in_set/enrich_U5U0_up$total_genes_in_set)*100
enrich_U5U0_up$perc_uni <- (enrich_U5U0_up$genes_in_universe/enrich_U5U0_up$total_genes_in_universe)*100
enrich_U5U0_down$perc_set <- (enrich_U5U0_down$genes_in_set/enrich_U5U0_down$total_genes_in_set)*-100
enrich_U5U0_down$perc_uni <- (enrich_U5U0_down$genes_in_universe/enrich_U5U0_down$total_genes_in_universe)*-100
enrich_U24U0_up$perc_set <- (enrich_U24U0_up$genes_in_set/enrich_U24U0_up$total_genes_in_set)*100
enrich_U24U0_up$perc_uni <- (enrich_U24U0_up$genes_in_universe/enrich_U24U0_up$total_genes_in_universe)*100
enrich_U24U0_down$perc_set <- (enrich_U24U0_down$genes_in_set/enrich_U24U0_down$total_genes_in_set)*-100
enrich_U24U0_down$perc_uni <- (enrich_U24U0_down$genes_in_universe/enrich_U24U0_down$total_genes_in_universe)*-100
enrich_T5T0_up$perc_set <- (enrich_T5T0_up$genes_in_set/enrich_T5T0_up$total_genes_in_set)*100
enrich_T5T0_up$perc_uni <- (enrich_T5T0_up$genes_in_universe/enrich_T5T0_up$total_genes_in_universe)*100
enrich_T5T0_down$perc_set <- (enrich_T5T0_down$genes_in_set/enrich_T5T0_down$total_genes_in_set)*-100
enrich_T5T0_down$perc_uni <- (enrich_T5T0_down$genes_in_universe/enrich_T5T0_down$total_genes_in_universe)*-100
enrich_T24T0_up$perc_set <- (enrich_T24T0_up$genes_in_set/enrich_T24T0_up$total_genes_in_set)*100
enrich_T24T0_up$perc_uni <- (enrich_T24T0_up$genes_in_universe/enrich_T24T0_up$total_genes_in_universe)*100
enrich_T24T0_down$perc_set <- (enrich_T24T0_down$genes_in_set/enrich_T24T0_down$total_genes_in_set)*-100
enrich_T24T0_down$perc_uni <- (enrich_T24T0_down$genes_in_universe/enrich_T24T0_down$total_genes_in_universe)*-100

enrich_U5U0_up$condition <- c('U5U0')
enrich_U5U0_down$condition <- c('U5U0')
enrich_U24U0_up$condition <- c('U24U0')
enrich_U24U0_down$condition <- c('U24U0')
enrich_T5T0_up$condition <- c('T5T0')
enrich_T5T0_down$condition <- c('T5T0')
enrich_T24T0_up$condition <- c('T24T0')
enrich_T24T0_down$condition <- c('T24T0')

enrich_U5U0_up$sign <- c('positive')
enrich_U5U0_down$sign <- c('negative')
enrich_U24U0_up$sign <- c('positive')
enrich_U24U0_down$sign <- c('negative')
enrich_T5T0_up$sign <- c('positive')
enrich_T5T0_down$sign <- c('negative')
enrich_T24T0_up$sign <- c('positive')
enrich_T24T0_down$sign <- c('negative')

enrich_U5U0_up_sig <- which(enrich_U5U0_up$adj_pvalue < 0.05)
enrich_U5U0_down_sig <- which(enrich_U5U0_down$adj_pvalue < 0.05)
enrich_U24U0_up_sig <- which(enrich_U24U0_up$adj_pvalue < 0.05)
enrich_U24U0_down_sig <- which(enrich_U24U0_down$adj_pvalue < 0.05)
enrich_T5T0_up_sig <- which(enrich_T5T0_up$adj_pvalue < 0.05)
enrich_T5T0_down_sig <- which(enrich_T5T0_down$adj_pvalue < 0.05)
enrich_T24T0_up_sig <- which(enrich_T24T0_up$adj_pvalue < 0.05)
enrich_T24T0_down_sig <- which(enrich_T24T0_down$adj_pvalue < 0.05)

enrich_sig <- unique(c(enrich_U5U0_up_sig,
                           enrich_U5U0_down_sig,
                           enrich_U24U0_up_sig,
                           enrich_U24U0_down_sig,
                           enrich_T5T0_up_sig,
                           enrich_T5T0_down_sig,
                           enrich_T24T0_up_sig,
                           enrich_T24T0_down_sig))

enrich_table <- rbind(enrich_U5U0_up[enrich_sig,],
                      enrich_U5U0_down[enrich_sig,],
                      enrich_U24U0_up[enrich_sig,],
                      enrich_U24U0_down[enrich_sig,],
                      enrich_T5T0_up[enrich_sig,],
                      enrich_T5T0_down[enrich_sig,],
                      enrich_T24T0_up[enrich_sig,],
                      enrich_T24T0_down[enrich_sig,])

enrich_table$cond_sign <- paste(enrich_table$condition,enrich_table$sign,sep='_')

enrich_table$term <- factor(enrich_table$term, levels = unique(sort(enrich_table$term, decreasing = TRUE)))

ggplot(enrich_table, aes(x = term)) +
  scale_y_continuous(limits = c(-25, 25), labels=abs) +
  geom_bar(aes(y = perc_set, fill = condition, color = condition), stat = "identity", position = position_dodge(0.75), width=0.5) +
  geom_bar(aes(y = perc_uni), color = "black", fill = NA, stat = "identity", position = "dodge",width=0.75) +
  #geom_text(data = enrich_table %>% filter(sign == "up"), aes(y = perc_set + 0.2, label = label), 
  #          size = 5) +
  #geom_text(data = enrich_table %>% filter(sign == "down"), aes(y = perc_set - 0.2, label = label), 
  #          size = 5) +
  coord_flip() + 
  scale_fill_manual(values = c("U5U0" = myPalette[7],
                               "U24U0" = myPalette[8],
                               "T5T0" = myPalette[4],
                               "T24T0" = myPalette[5]),
                    labels = c("U5","U24","P5","P24"),
                    limits = c("U5U0","U24U0","T5T0","T24T0")) +
  scale_color_manual(values = c("U5U0" = myPalette[7],
                                "U24U0" = myPalette[8],
                                "T5T0" = myPalette[4],
                                "T24T0" = myPalette[5]), guide = FALSE) +
  labs(x = "", y = "% of genes") +
  theme_bw() +
  theme(panel.grid.major.x = element_line("#CBCBCB", size = 0.5, linetype = "dotted"),
        text = element_text(size = 22))

##### heatmap for genes in the model ##### 

# find model directory
setwd('C:/Users/jmfic/OneDrive - University of Virginia/21_SUMMA/persister')
model_dir <- getwd()
model_dir <- paste0(model_dir, '/', 'model')

metabGenes <- read.csv(file.path(model_dir, 'data', 'iPau21genes.csv'))

# filter all DE data with for metab genes
T5T0_metab <- T5T0[which(T5T0$X %in% metabGenes$PA14_metabolic_model_genes),]
T24T0_metab <- T24T0[which(T24T0$X %in% metabGenes$PA14_metabolic_model_genes),]
U5U0_metab <- U5U0[which(U5U0$X %in% metabGenes$PA14_metabolic_model_genes),]
U24U0_metab <- U24U0[which(U24U0$X %in% metabGenes$PA14_metabolic_model_genes),]

T5T0_metab_sig <- filter(T5T0_metab,padj < 0.01)
U5U0_metab_sig <- filter(U5U0_metab,padj < 0.01)
T24T0_metab_sig <- filter(T24T0_metab,padj < 0.01)
U24U0_metab_sig <- filter(U24U0_metab,padj < 0.01)

metab_sig_genes<- unique(select(rbind(T5T0_metab_sig,U5U0_metab_sig,T24T0_metab_sig,U24U0_metab_sig),X))

heatmap_metab_sig <- metab_sig_genes
heatmap_metab_sig$X <- sort(heatmap_metab_sig$X)
heatmap_metab_sig$T5T0_metab <- rep(0,nrow(heatmap_metab_sig))
heatmap_metab_sig$U5U0_metab <- rep(0,nrow(heatmap_metab_sig))
heatmap_metab_sig$T24T0_metab <- rep(0,nrow(heatmap_metab_sig))
heatmap_metab_sig$U24U0_metab <- rep(0,nrow(heatmap_metab_sig))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
heatmap_metab_sig[heatmap_metab_sig$X %in% T5T0_metab$X,]$T5T0_metab = T5T0_metab[T5T0_metab$X %in% heatmap_metab_sig$X,]$log2FoldChange
heatmap_metab_sig[heatmap_metab_sig$X %in% U5U0_metab$X,]$U5U0_metab = U5U0_metab[U5U0_metab$X %in% heatmap_metab_sig$X,]$log2FoldChange
heatmap_metab_sig[heatmap_metab_sig$X %in% T24T0_metab$X,]$T24T0_metab = T24T0_metab[T24T0_metab$X %in% heatmap_metab_sig$X,]$log2FoldChange
heatmap_metab_sig[heatmap_metab_sig$X %in% U24U0_metab$X,]$U24U0_metab = U24U0_metab[U24U0_metab$X %in% heatmap_metab_sig$X,]$log2FoldChange

rownames(heatmap_metab_sig) <- heatmap_metab_sig$X
heatmap_metab_sig <- as.matrix(select(heatmap_metab_sig,-X))

myBreaks <- c(min(heatmap_metab_sig),seq(-2+1e-5,2-1e-5,length=99),max(heatmap_metab_sig))
hmcols <- c(colorRampPalette(c(myPalette[9],"white",myPalette[10]))(100))

colnames(heatmap_metab_sig) <- c('P5','U5','P24','U24')

pheatmap(heatmap_metab_sig,
         color = hmcols,
         breaks = myBreaks,
         show_rownames = FALSE,
         treeheight_row = 0,
         cellwidth = 25,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         fontsize = 20)

##### functional categories for genes in the model #####

T5T0_metab_sig_down <- filter(T5T0_metab_sig, log2FoldChange < 0)
T5T0_metab_sig_up <- filter(T5T0_metab_sig, log2FoldChange > 0)
U5U0_metab_sig_down <- filter(U5U0_metab_sig, log2FoldChange < 0)
U5U0_metab_sig_up <- filter(U5U0_metab_sig, log2FoldChange > 0)
T24T0_metab_sig_down <- filter(T24T0_metab_sig, log2FoldChange < 0)
T24T0_metab_sig_up <- filter(T24T0_metab_sig, log2FoldChange > 0)
U24U0_metab_sig_down <- filter(U24U0_metab_sig, log2FoldChange < 0)
U24U0_metab_sig_up <- filter(U24U0_metab_sig, log2FoldChange > 0)

metab_sig_down <- rbind(T5T0_metab_sig_down,
                        U5U0_metab_sig_down,
                        T24T0_metab_sig_down,
                        U24U0_metab_sig_down)

metab_sig_up <- rbind(T5T0_metab_sig_up,
                        U5U0_metab_sig_up,
                        T24T0_metab_sig_up,
                        U24U0_metab_sig_up)

# write.csv(metab_sig_down, file.path(dir, 'analysis', 'DESeq2', 'metab_sig_down.csv'))
# write.csv(metab_sig_up, file.path(dir, 'analysis', 'DESeq1', 'metab_sig_up.csv'))

##### find the expression profiles for a specific gene: #####
check <- filter(complete, X == "PA14_52820")

##### find genes that are significant in treated conditions but not untreated #####

# commonly downregulated:
T5T0_sig_down <- filter(T5T0,padj < 0.01 & log2FoldChange < 0)
U5U0_sig_down <- filter(U5U0,padj < 0.01 & log2FoldChange < 0)
T24T0_sig_down <- filter(T24T0,padj < 0.01 & log2FoldChange < 0)
U24U0_sig_down <- filter(U24U0,padj < 0.01 & log2FoldChange < 0)

treated_sig_down <- intersect(T5T0_sig_down$X,T24T0_sig_down$X)

treated_sig_down_unique <- setdiff(treated_sig_down, U5U0_sig_down$X)
treated_sig_down_unique <- setdiff(treated_sig_down_unique, U24U0_sig_down$X)

# commonly upregulated:
T5T0_sig_up <- filter(T5T0,padj < 0.01 & log2FoldChange > 0)
U5U0_sig_up <- filter(U5U0,padj < 0.01 & log2FoldChange > 0)
T24T0_sig_up <- filter(T24T0,padj < 0.01 & log2FoldChange > 0)
U24U0_sig_up <- filter(U24U0,padj < 0.01 & log2FoldChange > 0)

treated_sig_up <- intersect(T5T0_sig_up$X,T24T0_sig_up$X)

treated_sig_up_unique <- setdiff(treated_sig_up, U5U0_sig_up$X)
treated_sig_up_unique <- setdiff(treated_sig_up_unique, U24U0_sig_up$X)

treated_unique <- c(treated_sig_down_unique, treated_sig_up_unique)

write.csv(treated_unique, file='treated_unique.csv')
##### heatmap for succinate genes ##### 

succinateGenes <- c('PA14_03430','PA14_52670','PA14_52840','PA14_52820','PA14_52810','PA14_68290','PA14_62880','PA14_62860','PA14_44060','PA14_44050','PA14_44030','PA14_44020','PA14_43950','PA14_43940','PA14_49130','PA14_01460','PA14_38660','PA14_38640','PA14_72340','PA14_01460','PA14_49130','PA14_05230','PA14_52870','PA14_52630','PA14_49380','PA14_38660','PA14_30050','PA14_23930')

# filter all DE data with for metab genes
T5T0_succinate <- T5T0[which(T5T0$X %in% succinateGenes),]
T24T0_succinate <- T24T0[which(T24T0$X %in% succinateGenes),]
U5U0_succinate <- U5U0[which(U5U0$X %in% succinateGenes),]
U24U0_succinate <- U24U0[which(U24U0$X %in% succinateGenes),]

T5T0_succinate_sig <- T5T0_succinate
U5U0_succinate_sig <- U5U0_succinate
T24T0_succinate_sig <- T24T0_succinate
U24U0_succinate_sig <- U24U0_succinate

succinate_sig_genes<- unique(select(rbind(T5T0_succinate_sig,U5U0_succinate_sig,T24T0_succinate_sig,U24U0_succinate_sig),X))

heatmap_succinate_sig <- succinate_sig_genes
heatmap_succinate_sig$X <- sort(heatmap_succinate_sig$X)
heatmap_succinate_sig$T5T0_succinate <- rep(0,nrow(heatmap_succinate_sig))
heatmap_succinate_sig$U5U0_succinate <- rep(0,nrow(heatmap_succinate_sig))
heatmap_succinate_sig$T24T0_succinate <- rep(0,nrow(heatmap_succinate_sig))
heatmap_succinate_sig$U24U0_succinate <- rep(0,nrow(heatmap_succinate_sig))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
heatmap_succinate_sig[heatmap_succinate_sig$X %in% T5T0_succinate$X,]$T5T0_succinate = T5T0_succinate[T5T0_succinate$X %in% heatmap_succinate_sig$X,]$log2FoldChange
heatmap_succinate_sig[heatmap_succinate_sig$X %in% U5U0_succinate$X,]$U5U0_succinate = U5U0_succinate[U5U0_succinate$X %in% heatmap_succinate_sig$X,]$log2FoldChange
heatmap_succinate_sig[heatmap_succinate_sig$X %in% T24T0_succinate$X,]$T24T0_succinate = T24T0_succinate[T24T0_succinate$X %in% heatmap_succinate_sig$X,]$log2FoldChange
heatmap_succinate_sig[heatmap_succinate_sig$X %in% U24U0_succinate$X,]$U24U0_succinate = U24U0_succinate[U24U0_succinate$X %in% heatmap_succinate_sig$X,]$log2FoldChange

rownames(heatmap_succinate_sig) <- heatmap_succinate_sig$X
heatmap_succinate_sig <- as.matrix(select(heatmap_succinate_sig,-X))

myBreaks <- c(-2,seq(-2+1e-5,2-1e-5,length=99),2)
hmcols <- c(colorRampPalette(c(myPalette[9],"white",myPalette[10]))(100))

colnames(heatmap_succinate_sig) <- c('P5 vs. P0','U5 vs. U0','P24 vs. P0','U24 vs. U0')

pheatmap(heatmap_succinate_sig,
         color = hmcols,
         breaks = myBreaks,
         show_rownames = TRUE,
         treeheight_row = 0,
         treeheight_col = 0,
         cellwidth = 10,
         cellheight = 10,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         fontsize = 8,
         cluster_rows = FALSE
)


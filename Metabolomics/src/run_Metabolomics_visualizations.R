#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### set-up working environment #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# set the working directory
setwd('..')
dir <- getwd()

source(file.path(dir, 'src', "multiplot.R"))

# load the datasets
df_BIT_values <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_transformed.csv"))
df_BIT_uniStats_0 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_0.csv"))
df_BIT_uniStats_5 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_5.csv"))
df_BIT_uniStats_24 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_24.csv"))
df_BIT_uniStats_untreated <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_untreated.csv"))
df_BIT_uniStats_persister <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_persister.csv"))
df_BIT_uniStats_dead <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_dead.csv"))

df_BIT_values <- dplyr::select(df_BIT_values, -X)
df_BIT_uniStats_0 <- dplyr::select(df_BIT_uniStats_0, -X)
df_BIT_uniStats_5 <- dplyr::select(df_BIT_uniStats_5, -X)
df_BIT_uniStats_24 <- dplyr::select(df_BIT_uniStats_24, -X)
df_BIT_uniStats_untreated <- dplyr::select(df_BIT_uniStats_untreated, -X)
df_BIT_uniStats_persister <- dplyr::select(df_BIT_uniStats_persister, -X)
df_BIT_uniStats_dead <- dplyr::select(df_BIT_uniStats_dead, -X)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining overall significant metabolites for conditions #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# untreated significant metabs
untreated_metabs_5 <- filter(df_BIT_uniStats_untreated, condition == "untreated_5-0" & p < 0.05)
untreated_metabs_24 <- filter(df_BIT_uniStats_untreated, condition == "untreated_24-0" & p < 0.05)
untreated_metabs <- unique(c(untreated_metabs_5$metab,untreated_metabs_24$metab))

# persister significant metabs
persister_metabs_5 <- filter(df_BIT_uniStats_persister, condition == "persister_5-0" & p < 0.05)
persister_metabs_24 <- filter(df_BIT_uniStats_persister, condition == "persister_24-0" & p < 0.05)
persister_metabs <- unique(c(persister_metabs_5$metab,persister_metabs_24$metab))

# dead significant metabs
dead_metabs_5 <- filter(df_BIT_uniStats_dead, condition == "dead_5-0" & p < 0.05)
dead_metabs_24 <- filter(df_BIT_uniStats_dead, condition == "dead_24-0" & p < 0.05)
dead_metabs <- unique(c(dead_metabs_5$metab,dead_metabs_24$metab))

# to get subsystem information
df_untreated_metabs <- df_BIT_values[match(untreated_metabs, df_BIT_values$BIOCHEMICAL),]
df_persister_metabs <- df_BIT_values[match(persister_metabs, df_BIT_values$BIOCHEMICAL),]
df_dead_metabs <- df_BIT_values[match(dead_metabs, df_BIT_values$BIOCHEMICAL),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmaps for common metabolites #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# determine common significant metabs across conditions
common_metabs <- intersect(intersect(untreated_metabs, persister_metabs),dead_metabs)

# to get subsystem information
df_common_metabs <- df_BIT_values[match(common_metabs, df_BIT_values$BIOCHEMICAL),]

# generate common plots
untreated_metabs_vals_common <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% common_metabs,]
untreated_metabs_vals_common <- filter(untreated_metabs_vals_common, dose == 0)
persister_metabs_vals_common <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% common_metabs,]
persister_metabs_vals_common <- filter(persister_metabs_vals_common, dose == 0.1)
dead_metabs_vals_common <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% common_metabs,]
dead_metabs_vals_common <- filter(dead_metabs_vals_common, dose == 10)

untreated_metabs_mean0_common <- rep(0,length(common_metabs))
untreated_metabs_mean5_common <- rep(0,length(common_metabs))
untreated_metabs_mean24_common <- rep(0,length(common_metabs))
persister_metabs_mean0_common <- rep(0,length(common_metabs))
persister_metabs_mean5_common <- rep(0,length(common_metabs))
persister_metabs_mean24_common <- rep(0,length(common_metabs))
dead_metabs_mean0_common <- rep(0,length(common_metabs))
dead_metabs_mean5_common <- rep(0,length(common_metabs))
dead_metabs_mean24_common <- rep(0,length(common_metabs))
count <- 1

for (metab in common_metabs){
  df_metab0 <- filter(untreated_metabs_vals_common, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(untreated_metabs_vals_common, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(untreated_metabs_vals_common, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  untreated_metabs_mean0_common[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  untreated_metabs_mean5_common[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  untreated_metabs_mean24_common[count] <- metab_mean24
  
  df_metab0 <- filter(persister_metabs_vals_common, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(persister_metabs_vals_common, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(persister_metabs_vals_common, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  persister_metabs_mean0_common[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  persister_metabs_mean5_common[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  persister_metabs_mean24_common[count] <- metab_mean24
  
  df_metab0 <- filter(dead_metabs_vals_common, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(dead_metabs_vals_common, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(dead_metabs_vals_common, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  dead_metabs_mean0_common[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  dead_metabs_mean5_common[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  dead_metabs_mean24_common[count] <- metab_mean24
  
  count <- count + 1
}

untreated_metabs_mean_common <- data.frame(T0 = untreated_metabs_mean0_common,
                                    T5 = untreated_metabs_mean5_common,
                                    T24 = untreated_metabs_mean24_common)
rownames(untreated_metabs_mean_common) <- common_metabs
persister_metabs_mean_common <- data.frame(T0 = persister_metabs_mean0_common,
                                           T5 = persister_metabs_mean5_common,
                                           T24 = persister_metabs_mean24_common)
rownames(persister_metabs_mean_common) <- common_metabs
dead_metabs_mean_common <- data.frame(T0 = dead_metabs_mean0_common,
                                           T5 = dead_metabs_mean5_common,
                                           T24 = dead_metabs_mean24_common)
rownames(dead_metabs_mean_common) <- common_metabs


breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))

plot_untreated_metabs_mean_common <- pheatmap(untreated_metabs_mean_common,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

plot_persister_metabs_mean_common <- pheatmap(persister_metabs_mean_common,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

plot_dead_metabs_mean_common <- pheatmap(dead_metabs_mean_common,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

# re-make heatmaps so that they are all in the same order based on untreated clustering

row_order <- rownames(untreated_metabs_mean_common)[plot_untreated_metabs_mean_common$tree_row$order]
persister_metabs_mean_common <- persister_metabs_mean_common[match(row_order,row.names(persister_metabs_mean_common)),]
dead_metabs_mean_common <- dead_metabs_mean_common[match(row_order,row.names(dead_metabs_mean_common)),]

plot_persister_metabs_mean_common <- pheatmap(persister_metabs_mean_common,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              cluster_rows = FALSE)

plot_dead_metabs_mean_common <- pheatmap(dead_metabs_mean_common,
                                         color = cols,
                                         breaks = breaksList,
                                         border_color = "black",
                                         #show_rownames = FALSE,
                                         treeheight_row = 0,
                                         scale = "row",
                                         cellwidth = 15,
                                         cellheight = 8,
                                         cluster_cols = FALSE,
                                         cluster_rows = FALSE)
# 
# # re-make heatmaps so that they are all in the same order based on alpabetical SUB_PATHWAY
# row_order <- df_common_metabs[order(df_common_metabs$SUPER_PATHWAY),]$BIOCHEMICAL
# untreated_metabs_mean_common <- untreated_metabs_mean_common[match(row_order,row.names(untreated_metabs_mean_common)),]
# persister_metabs_mean_common <- persister_metabs_mean_common[match(row_order,row.names(persister_metabs_mean_common)),]
# dead_metabs_mean_common <- dead_metabs_mean_common[match(row_order,row.names(dead_metabs_mean_common)),]
# 
# plot_untreated_metabs_mean_common <- pheatmap(untreated_metabs_mean_common,
#                                               color = cols,
#                                               breaks = breaksList,
#                                               border_color = "black",
#                                               #show_rownames = FALSE,
#                                               treeheight_row = 0,
#                                               scale = "row",
#                                               cellwidth = 15,
#                                               cellheight = 8,
#                                               cluster_cols = FALSE,
#                                               cluster_rows = FALSE)
# 
# plot_persister_metabs_mean_common <- pheatmap(persister_metabs_mean_common,
#                                               color = cols,
#                                               breaks = breaksList,
#                                               border_color = "black",
#                                               #show_rownames = FALSE,
#                                               treeheight_row = 0,
#                                               scale = "row",
#                                               cellwidth = 15,
#                                               cellheight = 8,
#                                               cluster_cols = FALSE,
#                                               cluster_rows = FALSE)
# 
# plot_dead_metabs_mean_common <- pheatmap(dead_metabs_mean_common,
#                                          color = cols,
#                                          breaks = breaksList,
#                                          border_color = "black",
#                                          #show_rownames = FALSE,
#                                          treeheight_row = 0,
#                                          scale = "row",
#                                          cellwidth = 15,
#                                          cellheight = 8,
#                                          cluster_cols = FALSE,
#                                          cluster_rows = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmaps for unique metabolites #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# determine unique untreated significant metabs across conditions
untreated_unique_metabs <- setdiff(setdiff(untreated_metabs,persister_metabs),dead_metabs)
df_untreated_unique_metabs <- df_BIT_values[match(untreated_unique_metabs, df_BIT_values$BIOCHEMICAL),]

untreated_metabs_vals_unique <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% untreated_unique_metabs,]
untreated_metabs_vals_unique <- filter(untreated_metabs_vals_unique, dose == 0)

untreated_metabs_mean0_unique <- rep(0,length(untreated_unique_metabs))
untreated_metabs_mean5_unique <- rep(0,length(untreated_unique_metabs))
untreated_metabs_mean24_unique <- rep(0,length(untreated_unique_metabs))

count <- 1

for (metab in untreated_unique_metabs){
  df_metab0 <- filter(untreated_metabs_vals_unique, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(untreated_metabs_vals_unique, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(untreated_metabs_vals_unique, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  untreated_metabs_mean0_unique[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  untreated_metabs_mean5_unique[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  untreated_metabs_mean24_unique[count] <- metab_mean24
  
  count <- count + 1
}

untreated_metabs_mean_unique <- data.frame(T0 = untreated_metabs_mean0_unique,
                                           T5 = untreated_metabs_mean5_unique,
                                           T24 = untreated_metabs_mean24_unique)
rownames(untreated_metabs_mean_unique) <- untreated_unique_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))

plot_untreated_metabs_mean_unique <- pheatmap(untreated_metabs_mean_unique,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

# determine unique persister significant metabs across conditions
persister_unique_metabs <- setdiff(setdiff(persister_metabs,untreated_metabs),dead_metabs)
df_persister_unique_metabs <- df_BIT_values[match(persister_unique_metabs, df_BIT_values$BIOCHEMICAL),]

persister_metabs_vals_unique <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% persister_unique_metabs,]
persister_metabs_vals_unique <- filter(persister_metabs_vals_unique, dose == 0.1)

persister_metabs_mean0_unique <- rep(0,length(persister_unique_metabs))
persister_metabs_mean5_unique <- rep(0,length(persister_unique_metabs))
persister_metabs_mean24_unique <- rep(0,length(persister_unique_metabs))

count <- 1

for (metab in persister_unique_metabs){
  df_metab0 <- filter(persister_metabs_vals_unique, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(persister_metabs_vals_unique, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(persister_metabs_vals_unique, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  persister_metabs_mean0_unique[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  persister_metabs_mean5_unique[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  persister_metabs_mean24_unique[count] <- metab_mean24
  
  count <- count + 1
}

persister_metabs_mean_unique <- data.frame(T0 = persister_metabs_mean0_unique,
                                           T5 = persister_metabs_mean5_unique,
                                           T24 = persister_metabs_mean24_unique)
rownames(persister_metabs_mean_unique) <- persister_unique_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))

plot_persister_metabs_mean_unique <- pheatmap(persister_metabs_mean_unique,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

# determine unique dead significant metabs across conditions
dead_unique_metabs <- setdiff(setdiff(dead_metabs,untreated_metabs),persister_metabs)
df_dead_unique_metabs <- df_BIT_values[match(dead_unique_metabs, df_BIT_values$BIOCHEMICAL),]

dead_metabs_vals_unique <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% dead_unique_metabs,]
dead_metabs_vals_unique <- filter(dead_metabs_vals_unique, dose == 10)

dead_metabs_mean0_unique <- rep(0,length(dead_unique_metabs))
dead_metabs_mean5_unique <- rep(0,length(dead_unique_metabs))
dead_metabs_mean24_unique <- rep(0,length(dead_unique_metabs))

count <- 1

for (metab in dead_unique_metabs){
  df_metab0 <- filter(dead_metabs_vals_unique, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(dead_metabs_vals_unique, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(dead_metabs_vals_unique, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  dead_metabs_mean0_unique[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  dead_metabs_mean5_unique[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  dead_metabs_mean24_unique[count] <- metab_mean24
  
  count <- count + 1
}

dead_metabs_mean_unique <- data.frame(T0 = dead_metabs_mean0_unique,
                                           T5 = dead_metabs_mean5_unique,
                                           T24 = dead_metabs_mean24_unique)
rownames(dead_metabs_mean_unique) <- dead_unique_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))

plot_dead_metabs_mean_unique <- pheatmap(dead_metabs_mean_unique,
                                              color = cols,
                                              breaks = breaksList,
                                              border_color = "black",
                                              #show_rownames = FALSE,
                                              treeheight_row = 0,
                                              scale = "row",
                                              cellwidth = 15,
                                              cellheight = 8,
                                              cluster_cols = FALSE,
                                              clustering_distance_rows = "correlation",
                                              clustering_method = "average")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmap for complete untreated #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## metabs unique to untreated
untreated_metabs_vals <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% untreated_metabs,]
untreated_metabs_vals <- filter(untreated_metabs_vals, dose == 0)

untreated_metabs_mean0 <- rep(0,length(untreated_metabs))
untreated_metabs_mean5 <- rep(0,length(untreated_metabs))
untreated_metabs_mean24 <- rep(0,length(untreated_metabs))
count <- 1

for (metab in untreated_metabs){
  df_metab0 <- filter(untreated_metabs_vals, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(untreated_metabs_vals, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(untreated_metabs_vals, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  untreated_metabs_mean0[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  untreated_metabs_mean5[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  untreated_metabs_mean24[count] <- metab_mean24
  
  count <- count + 1
}

untreated_metabs_mean <- data.frame(T0 = untreated_metabs_mean0,
                                    T5 = untreated_metabs_mean5,
                                    T24 = untreated_metabs_mean24)
rownames(untreated_metabs_mean) <- untreated_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))
pheatmap(untreated_metabs_mean,
         color = cols,
         breaks = breaksList,
         border_color = "black",
         #show_rownames = FALSE,
         treeheight_row = 0,
         scale = "row",
         cellwidth = 7,
         cellheight = 2,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "average")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmap for complete peresister #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

persister_metabs_vals <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% persister_metabs,]
persister_metabs_vals <- filter(persister_metabs_vals, dose == 0.1)

persister_metabs_mean0 <- rep(0,length(persister_metabs))
persister_metabs_mean5 <- rep(0,length(persister_metabs))
persister_metabs_mean24 <- rep(0,length(persister_metabs))
count <- 1

for (metab in persister_metabs){
  df_metab0 <- filter(persister_metabs_vals, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(persister_metabs_vals, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(persister_metabs_vals, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  persister_metabs_mean0[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  persister_metabs_mean5[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  persister_metabs_mean24[count] <- metab_mean24
  
  count <- count + 1
}

persister_metabs_mean <- data.frame(T0 = persister_metabs_mean0,
                                    T5 = persister_metabs_mean5,
                                    T24 = persister_metabs_mean24)
rownames(persister_metabs_mean) <- persister_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))
pheatmap(persister_metabs_mean,
         color = cols,
         breaks = breaksList,
         border_color = "black",
         #show_rownames = FALSE,
         treeheight_row = 0,
         scale = "row",
         cellwidth = 7,
         cellheight = 7,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "average")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmap for complete dead #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dead_metabs_vals <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% dead_metabs,]
dead_metabs_vals <- filter(dead_metabs_vals, dose == 10)

dead_metabs_mean0 <- rep(0,length(dead_metabs))
dead_metabs_mean5 <- rep(0,length(dead_metabs))
dead_metabs_mean24 <- rep(0,length(dead_metabs))
count <- 1

for (metab in dead_metabs){
  df_metab0 <- filter(dead_metabs_vals, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(dead_metabs_vals, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(dead_metabs_vals, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  dead_metabs_mean0[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  dead_metabs_mean5[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  dead_metabs_mean24[count] <- metab_mean24
  
  count <- count + 1
}

dead_metabs_mean <- data.frame(T0 = dead_metabs_mean0,
                                    T5 = dead_metabs_mean5,
                                    T24 = dead_metabs_mean24)
rownames(dead_metabs_mean) <- dead_metabs

breaksList = seq(-1.5,1.5,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))
pheatmap(dead_metabs_mean,
         color = cols,
         breaks = breaksList,
         border_color = "black",
         #show_rownames = FALSE,
         treeheight_row = 0,
         scale = "row",
         cellwidth = 25,
         cellheight = 10,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "average")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### heatmap for superset #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

union_metabs <- union(union(untreated_metabs,persister_metabs),dead_metabs)

union_metabs_vals <- df_BIT_values[df_BIT_values$BIOCHEMICAL %in% union_metabs,]

untreated_metabs_vals_union <- filter(union_metabs_vals, dose == 0)
persister_metabs_vals_union <- filter(union_metabs_vals, dose == 0.1)
dead_metabs_vals_union <- filter(union_metabs_vals, dose == 10)

untreated_metabs_mean0_union <- rep(0,length(union_metabs))
untreated_metabs_mean5_union <- rep(0,length(union_metabs))
untreated_metabs_mean24_union <- rep(0,length(union_metabs))
persister_metabs_mean0_union <- rep(0,length(union_metabs))
persister_metabs_mean5_union <- rep(0,length(union_metabs))
persister_metabs_mean24_union <- rep(0,length(union_metabs))
dead_metabs_mean0_union <- rep(0,length(union_metabs))
dead_metabs_mean5_union <- rep(0,length(union_metabs))
dead_metabs_mean24_union <- rep(0,length(union_metabs))
count <- 1

for (metab in union_metabs){
  df_metab0 <- filter(untreated_metabs_vals_union, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(untreated_metabs_vals_union, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(untreated_metabs_vals_union, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  untreated_metabs_mean0_union[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  untreated_metabs_mean5_union[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  untreated_metabs_mean24_union[count] <- metab_mean24
  
  df_metab0 <- filter(persister_metabs_vals_union, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(persister_metabs_vals_union, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(persister_metabs_vals_union, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  persister_metabs_mean0_union[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  persister_metabs_mean5_union[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  persister_metabs_mean24_union[count] <- metab_mean24
  
  df_metab0 <- filter(dead_metabs_vals_union, BIOCHEMICAL == metab & time == 0)
  df_metab5 <- filter(dead_metabs_vals_union, BIOCHEMICAL == metab & time == 5)
  df_metab24 <- filter(dead_metabs_vals_union, BIOCHEMICAL == metab & time == 24)
  
  metab_mean0 <- mean(df_metab0$intensity)
  dead_metabs_mean0_union[count] <- metab_mean0
  
  metab_mean5 <- mean(df_metab5$intensity)
  dead_metabs_mean5_union[count] <- metab_mean5
  
  metab_mean24 <- mean(df_metab24$intensity)
  dead_metabs_mean24_union[count] <- metab_mean24
  
  count <- count + 1
}

metabs_mean_union <- data.frame(untreated_T0 = untreated_metabs_mean0_union,
                                untreated_T5 = untreated_metabs_mean5_union,
                                untreated_T24 = untreated_metabs_mean24_union,
                                persister_T0 = persister_metabs_mean0_union,
                                persister_T5 = persister_metabs_mean5_union,
                                persister_T24 = persister_metabs_mean24_union,
                                dead_T0 = dead_metabs_mean0_union,
                                dead_T5 = dead_metabs_mean5_union,
                                dead_T24 = dead_metabs_mean24_union)
rownames(metabs_mean_union) <- union_metabs


breaksList = seq(-3,3,by = 0.01)
cols <- colorRampPalette((brewer.pal(n = 7, name = "Purples")))(length(breaksList))

plot_metabs_mean_union <- pheatmap(metabs_mean_union,
                                  color = cols,
                                  breaks = breaksList,
                                  border_color = "black",
                                  #treeheight_row = 0,
                                  scale = "row",
                                  #cellwidth = 25,
                                  #cellheight = 10,
                                  #cluster_cols = FALSE,
                                  show_rownames = FALSE,
                                  clustering_distance_rows = "correlation",
                                  clustering_method = "average")


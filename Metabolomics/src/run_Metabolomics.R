#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### set-up environment ##### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(gplots)
library(plotly)
library(car)
library(userfriendlyscience)
library(nlme)
library(multcomp)

# set directory
setwd('C:/Users/jmf7ak/OneDrive - University of Virginia/21_SUMMA/persister/Metabolomics')
dir <- getwd()

# load functions
source(file.path(dir, 'src', 'multiplot.R'))

# load the dataset
df <- read.csv(file.path(dir, 'data', 'METABOLON_DATA_origscaled_R.csv'))

# tidy data
df_tidy <- df %>%
  select(-SUPER_PATHWAY, -SUB_PATHWAY, -COMP_ID, -PLATFORM, -CHEMICAL_ID, -RI, -MASS, -CAS, -PUBCHEM, -KEGG, -HMDB) %>% 
  gather(key="sample", value="intensity", 3:74) %>%
  # the following line parses the variable "sample" into 5 variables using "_" as the deliminator
  separate(sample, into=c("date","drug","time","dose","bioRep","techRep"), sep = "_")

df_tidy$bioRep <- paste(df_tidy$date,df_tidy$bioRep,sep="")

# # changes spaces in metabolite names to underscores
# df_tidy$BIOCHEMICAL <- sub(" ","_",df_tidy$BIOCHEMICAL)
# # change colon in metabolite names to Cs
# df_tidy$BIOCHEMICAL <- sub(":","C",df_tidy$BIOCHEMICAL)
# remove comma from intensity
df_tidy$intensity <- as.numeric(gsub(",","",df_tidy$intensity))

# get BIT dataset
df_BIT <- filter(df_tidy, drug == "BIT")

# get pathway information
# df_BIT$SUPER_PATHWAYS <- rep(df$SUPER_PATHWAY,72)
# df_BIT$SUB_PATHWAYS <- rep(df$SUB_PATHWAY,72)

# create a color palette
myPalette <- c("#000000","#999999", "#00c590", "#008a65", "#00503A", "#ffcd5d", "#E69F00", "#996900", "#0072B2", "#D55E00","#56B4E9", "#e2b2cc","#CC79A7","#b34481")

# black, gray, bright green, *green*[4], dark green, light yellow, *yellow*[7], dark yellow, blue, orange, light blue, *pink* [13], dark pink
# from here: http://www.colorhexa.com/e79f00
# and here: http://www.colorhexa.com/009e73
# and here: http://www.colorhexa.com/cc79a7

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### visualization of raw data #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rawIntensityViz <- function(dataframe, Metabolite) {
  
  df_0_0 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 0)
  p_0_0 <- ggplot(df_0_0, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0_0") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_0_5 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 5)
  p_0_5 <- ggplot(df_0_5, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0_5") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_0_24 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 24)
  p_0_24 <- ggplot(df_0_24, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0_24") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_0.1_0 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 0)
  p_0.1_0 <- ggplot(df_0.1_0, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0.1_5") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_0.1_5 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 5)
  p_0.1_5 <- ggplot(df_0.1_5, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0.1_5") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_0.1_24 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 24)
  p_0.1_24 <- ggplot(df_0.1_24, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("0.1_24") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_10_0 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 0)
  p_10_0 <- ggplot(df_10_0, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("10_5") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_10_5 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 5)
  p_10_5 <- ggplot(df_10_5, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("10_5") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  df_10_24 <- filter(dataframe, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 24)
  p_10_24 <- ggplot(df_10_24, aes(x = techRep, y = intensity, fill = bioRep)) +
    geom_bar(stat= "identity", position = "dodge") +
    ggtitle("10_24") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  #quartz()
  layout <- matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, byrow = TRUE)
  
  png(filename = file.path(dir, 'analysis', Metabolite), width = 800, height = 600)
  multiplot(p_0_0, p_0.1_0, p_10_0, p_0_5, p_0.1_5, p_10_5, p_0_24, p_0.1_24, p_10_24, layout = layout)
  dev.off()
}

rawIntensityViz(df_BIT, as.character(df_BIT$PATHWAY_SORTORDER[6]))

metabs <- unique(df_BIT$PATHWAY_SORTORDER)

for (i in 1:length(metabs)) {
  rawIntensityViz(df_BIT, as.character(df_BIT$PATHWAY_SORTORDER[i]))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### visualization of un-normalized averages #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

averageIntensityViz <- function(dateframe, Metabolite) {
  
  mean_0_0 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 0)$intensity)
  sd_0_0 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 0)$intensity)
  
  mean_0_5 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 5)$intensity)
  sd_0_5 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 5)$intensity)
  
  mean_0_24 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 24)$intensity)
  sd_0_24 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0, time == 24)$intensity)
  
  mean_0.1_0 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 0)$intensity)
  sd_0.1_0 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 0)$intensity)
  
  mean_0.1_5 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 5)$intensity)
  sd_0.1_5 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 5)$intensity)
  
  mean_0.1_24 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 24)$intensity)
  sd_0.1_24 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 0.1, time == 24)$intensity)
  
  mean_10_0 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 0)$intensity)
  sd_10_0 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 0)$intensity)
  
  mean_10_5 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 5)$intensity)
  sd_10_5 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 5)$intensity)
  
  mean_10_24 <- mean(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 24)$intensity)
  sd_10_24 <- sd(filter(df_BIT, PATHWAY_SORTORDER == Metabolite, dose == 10, time == 24)$intensity)
  
  df_mean <- data.frame(mean = c(mean_0_0, mean_0_5, mean_0_24, mean_0.1_0, mean_0.1_5, mean_0.1_24, mean_10_0, mean_10_5, mean_10_24),
                        stdev = c(sd_0_0, sd_0_5, sd_0_24, sd_0.1_0, sd_0.1_5, sd_0.1_24, sd_10_0, sd_10_5, sd_10_24),
                        condition = c("0_0", "0_5","0_24", "0.1_0", "0.1_5","0.1_24","10_0", "10_5","10_24"))
  
  levels(df_mean$condition) <- c("0_0", "0_5","0_24", "0.1_0", "0.1_5","0.1_24","10_0", "10_5","10_24")
  
  ggplot(df_mean, aes(x = condition, y = mean, color = condition)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean-stdev, ymax = mean+stdev), width = 0.2) +
    scale_color_manual(values = c(myPalette[7],myPalette[7],myPalette[7],myPalette[4],myPalette[4],myPalette[4],myPalette[13],myPalette[13],myPalette[13])) +
    ggtitle(Metabolite) + 
    ylab("Mean Intensity") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
}

averageIntensityViz(df_BIT, as.character(df_BIT$PATHWAY_SORTORDER[6]))

metabs <- unique(df_BIT$PATHWAY_SORTORDER)

for (i in 1:length(metabs)) {
  averageIntensityViz(df_BIT, as.character(df_BIT$PATHWAY_SORTORDER[i]))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PCA of raw data (filtered out metabolites containing NAs) #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to remove all NAs...
df_BIT_NA <- within(df_BIT,  id <- paste(dose, time, bioRep, techRep, sep="_"))
df_BIT_NA <- select(df_BIT_NA, BIOCHEMICAL, intensity, id)
df_BIT_NA <- spread(df_BIT_NA, id, intensity)
rownames(df_BIT_NA) <- df_BIT_NA[,1]
df_BIT_NA[,1] <- NULL

df_BIT_NA <- drop_na(df_BIT_NA)

df_BIT_NA_pca <- t(df_BIT_NA)

# perform PCA
NA_pca_res <- prcomp(df_BIT_NA_pca)

plot(NA_pca_res)
summary(NA_pca_res)

NA_scores_BIT <- as.data.frame(NA_pca_res$x)

g <- ggplot(data = NA_scores_BIT, aes(x = PC1, y = PC2, label = rownames(NA_scores_BIT), fill = rownames(NA_scores_BIT))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour = c(myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],
                        myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],
                        myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],
                        myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],
                        myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],
                        myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],
                        myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],
                        myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],
                        myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14]), 
             alpha = 1, 
             size = 4) +
  scale_y_continuous(name = paste("PC2",":",as.character(summary(NA_pca_res)$importance[2,][2]*100),"%")) +
  scale_x_continuous(name = paste("PC1",":",as.character(summary(NA_pca_res)$importance[2,][1]*100),"%")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = 'none')

g
ggplotly(g)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Data imputation #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to impute with half the minimum...
# # should perhaps consider discarding the metabolites for which I need to impute values for an entire condition (i.e., if none of the time 0 samples detected a certain metabolite, should i still consider it?)
df_BIT_impute <- df_BIT

metabs <- unique(df_BIT$PATHWAY_SORTORDER)

missing_metabs <- c()

for (i in 1:length(metabs)) {
  
  df_metabolite <- filter(df_BIT_impute, PATHWAY_SORTORDER == as.character(df_BIT_impute$PATHWAY_SORTORDER[i]))
  
  if (all(is.na(filter(df_metabolite, time == 0 & dose == 0)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 5 & dose == 0)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 24 & dose == 0)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 0 & dose == 0.1)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 5 & dose == 0.1)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 24 & dose == 0.1)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 0 & dose == 10)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 5 & dose == 10)$intensity)) == FALSE &
      all(is.na(filter(df_metabolite, time == 24 & dose == 10)$intensity)) == FALSE) {
    minimum <- as.numeric(min(df_metabolite$intensity, na.rm=TRUE))
    df_metabolite$intensity[is.na(df_metabolite$intensity)] <- minimum/2
    df_BIT_impute[which(df_BIT_impute$PATHWAY_SORTORDER == as.character(df_BIT_impute$PATHWAY_SORTORDER[i])),] <- df_metabolite
  } else{
    missing_metabs <- append(missing_metabs, i)
  }
}

for (i in 1:length(missing_metabs)) {
  if (i == 1){
    df_BIT_impute <- subset(df_BIT_impute, PATHWAY_SORTORDER!=df_BIT_impute$PATHWAY_SORTORDER[missing_metabs[i]])
  } else {
    df_BIT_impute <- subset(df_BIT_impute, PATHWAY_SORTORDER!=df_BIT_impute$PATHWAY_SORTORDER[missing_metabs[i]-(i-1)])
  }
}

# # to save the imputed df
# df_BIT_impute_save <- df_BIT_impute
# df_BIT_impute_save <- within(df_BIT_impute_save,  id <- paste(dose, time, bioRep, techRep, sep="_"))
# df_BIT_impute_save$SUPER_PATHWAY <- rep(df$SUPER_PATHWAY[-missing_metabs],72)
# df_BIT_impute_save$SUB_PATHWAY <- rep(df$SUB_PATHWAY[-missing_metabs],72)
# write.csv(df_BIT_impute_save,'METABOLON_DATA_imputed.csv')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PCA of imputed, un-normalized samples #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_impute_pca <- within(df_BIT_impute,  id <- paste(dose, time, bioRep, techRep, sep="_"))
df_BIT_impute_pca <- select(df_BIT_impute_pca, BIOCHEMICAL, intensity, id)
df_BIT_impute_pca <- spread(df_BIT_impute_pca, id, intensity)
rownames(df_BIT_impute_pca) <- df_BIT_impute_pca[,1]
df_BIT_impute_pca[,1] <- NULL

df_BIT_impute_pca <- t(df_BIT_impute_pca)

# perform PCA
impute_pca_res <- prcomp(df_BIT_impute_pca)

plot(impute_pca_res)
summary(impute_pca_res)

impute_scores_BIT <- as.data.frame(impute_pca_res$x)

g <- ggplot(data = impute_scores_BIT, aes(x = PC1, y = PC2, label = rownames(impute_scores_BIT), fill = rownames(impute_scores_BIT))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour = c(myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],
                        myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],
                        myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],
                        myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],
                        myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],
                        myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],
                        myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],
                        myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],
                        myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14]), 
            alpha = 1, 
            size = 4) +
  scale_y_continuous(name = paste("PC2",":",as.character(summary(impute_pca_res)$importance[2,][2]*100),"%")) +
  scale_x_continuous(name = paste("PC1",":",as.character(summary(impute_pca_res)$importance[2,][1]*100),"%")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = 'none')

g
ggplotly(g)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Data transformation #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Reference: Wanichthanarak, K. et al. 2017. Metabox: A toolbox for metabolomic data analysis, interpretation and integrative exploration
# this reference discusses different aspects of processing metabolomics data to get it ready for univariate and multivariate analyses
# I particularly use this reference for the sample normalization

## Reference: van den Berg, R. et al. 2006. Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC Genomics
# this reference discusses different aspects of processing metabolomics data to get it ready for univariate and multivariate analyses

# data is in need of processing as evidenced by this skewed dataset
hist(df_BIT_impute$intensity, breaks = 100)

df_BIT_trans <- df_BIT_impute

## Modify the dataframe so that is it structured appropriately
# collapse the sample info to variable ID
df_BIT_trans <- within(df_BIT_trans,  id <- paste(dose, time, bioRep, techRep, sep="_"))
# now that ID is created, remove all the other sample info columns
df_BIT_trans <- select(df_BIT_trans, BIOCHEMICAL, intensity, id)
# make the dataframe wide
df_BIT_trans <- spread(df_BIT_trans, id, intensity)
# make the rownames the metabolite IDs
rownames(df_BIT_trans) <- df_BIT_trans[,1]
# get rid of the metabolite ID column
df_BIT_trans[,1] <- NULL

## Sample normalization
# Within each sample, I am dividing each metabolite peak intensity by the sum of all the intensities for that sample
# This step may not be necessary...it is unclear from UC Davis
# conditions <- colnames(df_BIT_trans)
# 
# for (i in 1:length(conditions)) {
#   df_BIT_trans[,i] <- df_BIT_trans[,i]/sum(df_BIT_trans[,i])
# }

# transpose the dataframe and ensure it is a dataframe and not a matrix
df_BIT_trans <- as.data.frame(t(df_BIT_trans))

## Transform the dataframe
# log2 transform because from a biological standpoint, it is easy to interpret doublings
# the log2 transform is commonly used across metabolomics experiments/papers in order to adjust the skewed data
df_BIT_trans <- log2(df_BIT_trans)
# visualize how log transforming impact the data, more normal but still skewed
hist(as.matrix(df_BIT_trans), breaks = 100)

## Center and auto-scale the dataset
# Centering converts all the fluctuations around zero instead of around the mean of the metabolite concentrations. Meant to focus on the fluctuating part of the data, leaving only the relevant variation
# Scaling aims to adjust for the differences in fold differences between the different metabolites. It helps to make all metabolites equally important
# Here, we apply autoscaling, which uses the standard deviation as the scaling factor. After autoscaling, all metabolites have a standard deviation of one and therefore the data is analyzed on the basis of correlations instead of covariances
df_BIT_trans <- scale(df_BIT_trans, center = TRUE, scale = FALSE)
# visualize how centering and scaling impacted the data, now more normally distributed
hist(df_BIT_trans, breaks = 100)

# # to save the transformed df
# df_BIT_trans_save <- as.data.frame(df_BIT_trans)
# df_BIT_trans_save$sample <- rownames(df_BIT_trans_save)
# rownames(df_BIT_trans_save) <- c()
# df_BIT_trans_save <- separate(df_BIT_trans_save, sample, into=c("dose","time","bioRep","techRep"), sep = "_")
# df_BIT_trans_save <- gather(df_BIT_trans_save,key="BIOCHEMICAL",value = "intensity",-dose,-time,-bioRep,-techRep)
# df_BIT_trans_save$SUPER_PATHWAY <- rep(df$SUPER_PATHWAY[-missing_metabs],72)
# df_BIT_trans_save$SUB_PATHWAY <- rep(df$SUB_PATHWAY[-missing_metabs],72)
# write.csv(df_BIT_save,'METABOLON_DATA_transformed.csv')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PCA of imputed, normalized samples #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# perform PCA
trans_pca_res <- prcomp(df_BIT_trans)

plot(trans_pca_res)
summary(trans_pca_res)

trans_scores_BIT <- as.data.frame(trans_pca_res$x)

g <- ggplot(data = trans_scores_BIT, aes(x = PC1, y = PC2, label = rownames(trans_scores_BIT), fill = rownames(trans_scores_BIT))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour = c(myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],myPalette[6],
                        myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],myPalette[7],
                        myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],myPalette[8],
                        myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],myPalette[3],
                        myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],myPalette[5],
                        myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],myPalette[4],
                        myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],myPalette[12],
                        myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],myPalette[14],
                        myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13],myPalette[13]), 
             alpha = 1, 
             size = 4) +
  scale_y_continuous(name = paste("PC2",":",as.character(summary(trans_pca_res)$importance[2,][2]*100),"%")) +
  scale_x_continuous(name = paste("PC1",":",as.character(summary(trans_pca_res)$importance[2,][1]*100),"%")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = 'none')

g
ggplotly(g)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### log2FC calculation on imputed data #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Reference: Vinaiza, M. et al. 2012. A guideline to univariate statistical analysis for LC/MS-based untargeted metabolomics-derived data. Metabolites.
# this reference indicates that raw intensity values should be used to calculate fold change and sets the standard of importance as 2 but admits this value is arbitrary

df_BIT_FC <- df_BIT_impute

## Modify the dataframe so that is it structured appropriately
# collapse the sample info to variable ID
df_BIT_FC <- within(df_BIT_FC,  id <- paste(dose, time, bioRep, techRep, sep="_"))
# now that ID is created, remove all the other sample info columns
df_BIT_FC <- select(df_BIT_FC, BIOCHEMICAL, intensity, id)
# make the dataframe wide
df_BIT_FC <- spread(df_BIT_FC, id, intensity)
# make the rownames the metabolite IDs
rownames(df_BIT_FC) <- df_BIT_FC[,1]
# get rid of the metabolite ID column
df_BIT_FC[,1] <- NULL
# transpose the dataframe and ensure it is a dataframe and not a matrix
df_BIT_FC <- as.data.frame(t(df_BIT_FC))
# create a new variable, sample, that is the ID names
df_BIT_FC$sample <- rownames(df_BIT_FC)
# remove row names
rownames(df_BIT_FC) <- c()
# separate the IDs for ease of calling
df_BIT_FC <- separate(df_BIT_FC, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_FC)[1:336]

# create a daframe to save the log2FC results
df_FC_sample <- data.frame(matrix(NA, nrow = 6, ncol = length(metabs)))
row.names(df_FC_sample) <- c("0vs0_5", "0vs0_24",
                             "0.1vs0.1_5","0.1vs0.1_24",
                             "10vs10_5","10vs10_24")
colnames(df_FC_sample) <- metabs

df_FC_time <- data.frame(matrix(NA, nrow = 6, ncol = length(metabs)))
row.names(df_FC_time) <- c("0v0.1_5","0v0.1_24",
                           "0v10_5","0v10_24",
                           "0.1v10_5","0.1v10_24")
colnames(df_FC_time) <- metabs

# calculate log2FC for sample constant conditions
for (i in 1:length(metabs)) {
  df_FC_sample[1,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5]))
  }
  df_FC_sample[2,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24]))
  }
  df_FC_sample[3,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5]))
  }
  df_FC_sample[4,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24]))
  }
  df_FC_sample[5,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]))
  }
  df_FC_sample[6,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 0])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]))
  }
}

df_FC_sample <- as.data.frame(t(df_FC_sample))

# calculate log2FC for time constant conditions
for (i in 1:length(metabs)) {
  df_FC_time[1,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5]))
  }
  df_FC_time[2,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24]))
  }
  df_FC_time[3,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]))
  }
  df_FC_time[4,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]))
  }
  df_FC_time[5,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 5])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 5]))
  }
  df_FC_time[6,i] <- if(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]) > mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24])){
    log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24]))
  }else{
    -1*log2(mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 0.1 & df_BIT_FC$time == 24])/mean(df_BIT_FC[,metabs[i]][df_BIT_FC$dose == 10 & df_BIT_FC$time == 24]))
  }
}

df_FC_time <- as.data.frame(t(df_FC_time))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### log2FC calculation on transformed data #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_trans_FC <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_trans_FC$sample <- rownames(df_BIT_trans_FC)
# remove row names
rownames(df_BIT_trans_FC) <- c()
# separate the IDs for ease of calling
df_BIT_trans_FC <- separate(df_BIT_trans_FC, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_trans_FC)[1:336]

# create a daframe to save the log2FC results
df_trans_FC_sample <- data.frame(matrix(NA, nrow = 9, ncol = length(metabs)))
row.names(df_trans_FC_sample) <- c("untreated_24-0", "untreated_5-0", "untreated_5-24",
                                   "persister_24-0", "persister_5-0", "persister_5-24",
                                   "dead_24-0", "dead_5-0", "dead_5-24")
colnames(df_trans_FC_sample) <- metabs

df_trans_FC_time <- data.frame(matrix(NA, nrow = 9, ncol = length(metabs)))
row.names(df_trans_FC_time) <- c("0_0.1-0","0_10-0","0_10-0.1",
                                 "5_0.1-0","5_10-0","5_10-0.1",
                                 "24_0.1-0","24_10-0","24_10-0.1")
colnames(df_trans_FC_time) <- metabs

# calculate log2FC for sample constant conditions
for (i in 1:length(metabs)) {
  
  ## untreated_24-0
  df_trans_FC_sample[1,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])))
  }
  
  ## untreated_5-0
  df_trans_FC_sample[2,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])))
  }
  
  ## untreated_5-24
  df_trans_FC_sample[3,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])))
  }
  
  ## persister_24-0
  df_trans_FC_sample[4,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])))
  }
  
  ## persister_5-0
  df_trans_FC_sample[5,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])))
  }
  
  ## persister_5-24
  df_trans_FC_sample[6,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])))
  }
  
  ## dead_24-0
  df_trans_FC_sample[7,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])))
  }
  
  ## dead_5-0
  df_trans_FC_sample[8,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])))
  }
  
  ## dead_5-24
  df_trans_FC_sample[9,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])))
  }
}

df_trans_FC_sample <- as.data.frame(t(df_trans_FC_sample))

# calculate log2FC for time constant conditions
for (i in 1:length(metabs)) {
  
  ## 0_0.1-0
  df_trans_FC_time[1,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])))
  }
  
  ## 0_10-0
  df_trans_FC_time[2,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])))
  }
  
  ## 0_10-0.1
  df_trans_FC_time[3,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 0])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 0])))
  }
  
  ## 5_0.1-0
  df_trans_FC_time[4,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])))
  }
  
  ## 5_10-0
  df_trans_FC_time[5,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])))
  }
  
  ## 5_10-0.1
  df_trans_FC_time[6,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 5])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 5])))
  }
  
  ## 24_0.1-0
  df_trans_FC_time[7,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])))
  }
  
  ## 24_10-0
  df_trans_FC_time[8,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])))
  }
  
  ## 24_10-0.1
  df_trans_FC_time[9,i] <- if(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24]) > mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])){
    log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])))
  }else{
    -1*log2(abs(mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 0.1 & df_BIT_trans_FC$time == 24])/mean(df_BIT_trans_FC[,metabs[i]][df_BIT_trans_FC$dose == 10 & df_BIT_trans_FC$time == 24])))
  }
}

df_trans_FC_time <- as.data.frame(t(df_trans_FC_time))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: univariate assumptions assessment #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## Normality
# because our sample size is small (8 replicates), statistical normality tests like Shapiro-Wilks are not really appropriate
# we can either assume normality (and run a one-way ANOVA) OR we can perform a non-parametric test (Kruskal-Wallis H Test)

## Homogeneity of variances
# use Levene's Test, which will compare the variances between two samples, so will need to set up different pairwise comparisons
# *MOST* of the samples have equal variances, which means that a more conservative approach would be better to use (e.g., Welch's ANOVA and Games-Howell post hoc test)

df_levene_0 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_levene_0) <- c("p-val")
colnames(df_levene_0) <- metabs

df_levene_5 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_levene_5) <- c("p-val")
colnames(df_levene_5) <- metabs

df_levene_24 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_levene_24) <- c("pval")
colnames(df_levene_24) <- metabs

for (i in 1:length(metabs)) {
  # time 0
  data_0 <- data.frame(intensity = df_BIT_uni[,metabs[i]][df_BIT_uni$time == 0],
                       dose = df_BIT_uni$dose[df_BIT_uni$time == 0])
  df_levene_0[i] <- leveneTest(intensity ~ dose, data = data_0)[3]
  
  # time 5
  data_5 <- data.frame(intensity = df_BIT_uni[,metabs[i]][df_BIT_uni$time == 5],
                       dose = df_BIT_uni$dose[df_BIT_uni$time == 5])
  df_levene_5[i] <- leveneTest(intensity ~ dose, data = data_5)[3]
  
  # time 24
  data_24 <- data.frame(intensity = df_BIT_uni[,metabs[i]][df_BIT_uni$time == 24],
                       dose = df_BIT_uni$dose[df_BIT_uni$time == 24])
  df_levene_24[i] <- leveneTest(intensity ~ dose, data = data_24)[3]
}

length(which(df_levene_0 < 0.05)) # 12
length(which(df_levene_5 < 0.05)) # 71
length(which(df_levene_24 < 0.05)) # 101

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: Welch one-way ANOVA analysis with Games-Howell post-hoc test #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## One-way Welch ANOVA for the different time-points
df_BIT_welch <- df_BIT_uni
df_BIT_welch <- select(df_BIT_welch, -techRep)

# create a daframe to save the time Welch ANOVA results
df_welch_0 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_welch_0) <- c("0")
colnames(df_welch_0) <- metabs

df_welch_5 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_welch_5) <- c("5")
colnames(df_welch_5) <- metabs

df_welch_24 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_welch_24) <- c("24")
colnames(df_welch_24) <- metabs

df_games_0 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_0) <- c("0_0.1-0","0_10-0","0_10-0.1")
colnames(df_games_0) <- metabs

df_games_5 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_5) <- c("5_0.1-0","5_10-0","5_10-0.1")
colnames(df_games_5) <- metabs

df_games_24 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_24) <- c("24_0.1-0","24_10-0","24_10-0.1")
colnames(df_games_24) <- metabs

# run Welch ANOVA for all the 0-hour, 5-hour, and 24-hour conditions
for (i in 1:length(metabs)) {
  
  ## time 0
  data_0 <- data.frame(intensity = df_BIT_welch[,metabs[i]][df_BIT_welch$time == 0],
                       dose = df_BIT_welch$dose[df_BIT_welch$time == 0])
  lm_0 <- oneway.test(intensity~dose, data = data_0, var.equal = FALSE)
  df_welch_0[1,i] <- lm_0$p.value
  posthoc_0 <- posthocTGH(data_0$intensity, data_0$dose, method=c("games-howell"), digits = 5)
  df_games_0[1,i] <- unlist(posthoc_0$output[2])[16]
  df_games_0[2,i] <- unlist(posthoc_0$output[2])[17]
  df_games_0[3,i] <- unlist(posthoc_0$output[2])[18]
  
  ## time 5
  data_5 <- data.frame(intensity = df_BIT_welch[,metabs[i]][df_BIT_welch$time == 5],
                       dose = df_BIT_welch$dose[df_BIT_welch$time == 5])
  lm_5 <- oneway.test(intensity~dose, data = data_5, var.equal = FALSE)
  df_welch_5[1,i] <- lm_5$p.value
  posthoc_5 <- posthocTGH(data_5$intensity, data_5$dose, method=c("games-howell"), digits = 5)
  df_games_5[1,i] <- unlist(posthoc_5$output[2])[16]
  df_games_5[2,i] <- unlist(posthoc_5$output[2])[17]
  df_games_5[3,i] <- unlist(posthoc_5$output[2])[18]
  
  ## time 24
  data_24 <- data.frame(intensity = df_BIT_welch[,metabs[i]][df_BIT_welch$time == 24],
                        dose = df_BIT_welch$dose[df_BIT_welch$time == 24])
  lm_24 <- oneway.test(intensity~dose, data = data_24, var.equal = FALSE)
  df_welch_24[1,i] <- lm_24$p.value
  posthoc_24 <- posthocTGH(data_24$intensity, data_24$dose, method=c("games-howell"), digits = 5)
  df_games_24[1,i] <- unlist(posthoc_24$output[2])[16]
  df_games_24[2,i] <- unlist(posthoc_24$output[2])[17]
  df_games_24[3,i] <- unlist(posthoc_24$output[2])[18]
}

# Multiple comparison testing using Benjamini Hochberg
# It is necessary to correct for multiple comparisons AFTER running ANOVA
df_welch_0_padj <- as.data.frame(t(df_welch_0))
df_welch_5_padj <- as.data.frame(t(df_welch_5))
df_welch_24_padj <- as.data.frame(t(df_welch_24))
df_welch_0_padj[,1] <- p.adjust(df_welch_0_padj[,1], method = "BH")
df_welch_5_padj[,1] <- p.adjust(df_welch_5_padj[,1], method = "BH")
df_welch_24_padj[,1] <- p.adjust(df_welch_24_padj[,1], method = "BH")

# combine anova and post-hoc results and filter for significant anova results
df_welch_res_0 <- cbind(t(df_games_0),df_welch_0_padj)
colnames(df_welch_res_0) <- c("0_0.1-0","0_10-0","0_10-0.1","p.adj")
df_welch_res_0 <- subset(df_welch_res_0, p.adj < 0.05)

df_welch_res_5 <- cbind(t(df_games_5),df_welch_5_padj)
colnames(df_welch_res_5) <- c("5_0.1-0","5_10-0","5_10-0.1","p.adj")
df_welch_res_5 <- subset(df_welch_res_5, p.adj < 0.05)

df_welch_res_24 <- cbind(t(df_games_24),df_welch_24_padj)
colnames(df_welch_res_24) <- c("24_0.1-0","24_10-0","24_10-0.1","p.adj")
df_welch_res_24 <- subset(df_welch_res_24, p.adj < 0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: merge welch anova-log2fc #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_FC_0 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_welch_res_0))
df_FC_5 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_welch_res_5))
df_FC_24 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_welch_res_24))

df_BIT_0_0.1v0 <- data.frame(df_FC_0[1], df_welch_res_0[1])
colnames(df_BIT_0_0.1v0) <- c("log2FC","p")
df_BIT_0_0.1v0$condition <- c("0_0.1-0")
df_BIT_0_0.1v0$metab <- rownames(df_BIT_0_0.1v0)

df_BIT_0_10v0 <- data.frame(df_FC_0[2], df_welch_res_0[2])
colnames(df_BIT_0_10v0) <- c("log2FC","p")
df_BIT_0_10v0$condition <- c("0_10-0")
df_BIT_0_10v0$metab <- rownames(df_BIT_0_10v0)

df_BIT_0_10v0.1 <- data.frame(df_FC_0[3], df_welch_res_0[3])
colnames(df_BIT_0_10v0.1) <- c("log2FC","p")
df_BIT_0_10v0.1$condition <- c("0_10-0.1")
df_BIT_0_10v0.1$metab <- rownames(df_BIT_0_10v0.1)

df_BIT_5_0.1v0 <- data.frame(df_FC_5[4], df_welch_res_5[1])
colnames(df_BIT_5_0.1v0) <- c("log2FC","p")
df_BIT_5_0.1v0$condition <- c("5_0.1-0")
df_BIT_5_0.1v0$metab <- rownames(df_BIT_5_0.1v0)

df_BIT_5_10v0 <- data.frame(df_FC_5[5], df_welch_res_5[2])
colnames(df_BIT_5_10v0) <- c("log2FC","p")
df_BIT_5_10v0$condition <- c("5_10-0")
df_BIT_5_10v0$metab <- rownames(df_BIT_5_10v0)

df_BIT_5_10v0.1 <- data.frame(df_FC_5[6], df_welch_res_5[3])
colnames(df_BIT_5_10v0.1) <- c("log2FC","p")
df_BIT_5_10v0.1$condition <- c("5_10-0.1")
df_BIT_5_10v0.1$metab <- rownames(df_BIT_5_10v0.1)

df_BIT_24_0.1v0 <- data.frame(df_FC_24[7], df_welch_res_24[1])
colnames(df_BIT_24_0.1v0) <- c("log2FC","p")
df_BIT_24_0.1v0$condition <- c("24_0.1-0")
df_BIT_24_0.1v0$metab <- rownames(df_BIT_24_0.1v0)

df_BIT_24_10v0 <- data.frame(df_FC_24[8], df_welch_res_24[2])
colnames(df_BIT_24_10v0) <- c("log2FC","p")
df_BIT_24_10v0$condition <- c("24_10-0")
df_BIT_24_10v0$metab <- rownames(df_BIT_24_10v0)

df_BIT_24_10v0.1 <- data.frame(df_FC_24[9], df_welch_res_24[3])
colnames(df_BIT_24_10v0.1) <- c("log2FC","p")
df_BIT_24_10v0.1$condition <- c("24_10-0.1")
df_BIT_24_10v0.1$metab <- rownames(df_BIT_24_10v0.1)

df_BIT_uniStats_0 <- rbind(df_BIT_0_0.1v0, df_BIT_0_10v0, df_BIT_0_10v0.1)
df_BIT_uniStats_0$condition <- factor(df_BIT_uniStats_0$condition,levels=c("0_0.1-0","0_10-0","0_10-0.1"))

df_BIT_uniStats_5 <- rbind(df_BIT_5_0.1v0, df_BIT_5_10v0, df_BIT_5_10v0.1)
df_BIT_uniStats_5$condition <- factor(df_BIT_uniStats_5$condition,levels=c("5_0.1-0","5_10-0","5_10-0.1"))

df_BIT_uniStats_24 <- rbind(df_BIT_24_0.1v0, df_BIT_24_10v0, df_BIT_24_10v0.1)
df_BIT_uniStats_24$condition <- factor(df_BIT_uniStats_24$condition,levels=c("24_0.1-0","24_10-0","24_10-0.1"))

write.csv(df_BIT_uniStats_0, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_0.csv'))
write.csv(df_BIT_uniStats_5, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_5.csv'))
write.csv(df_BIT_uniStats_24, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_24.csv'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: one-way ANOVA analysis with Tukey post-hoc test #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## One-way ANOVA for the different time-points
df_BIT_anova <- df_BIT_uni
df_BIT_anova <- select(df_BIT_anova, -techRep)

# create a daframe to save the time ANOVA results
df_anova_0 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_0) <- c("0")
colnames(df_anova_0) <- metabs

df_anova_5 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_5) <- c("5")
colnames(df_anova_5) <- metabs

df_anova_24 <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_24) <- c("24")
colnames(df_anova_24) <- metabs

df_tukey_0 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_0) <- c("0_0.1-0","0_10-0","0_10-0.1")
colnames(df_tukey_0) <- metabs

df_tukey_5 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_5) <- c("5_0.1-0","5_10-0","5_10-0.1")
colnames(df_tukey_5) <- metabs

df_tukey_24 <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_24) <- c("24_0.1-0","24_10-0","24_10-0.1")
colnames(df_tukey_24) <- metabs

# run ANOVA for all the 0-hour, 5-hour, and 24-hour conditions
for (i in 1:length(metabs)) {
  
  ## time 0
  data_0 <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$time == 0],
                       dose = df_BIT_anova$dose[df_BIT_anova$time == 0])
  lm_0 <- lm(formula = intensity ~ dose, data = data_0)
  df_anova_0[1,i] <- anova(lm_0)$`Pr(>F)`[1]
  a_0 <- aov(data_0$intensity ~ data_0$dose)
  posthoc_0 <- TukeyHSD(x = a_0, 'data_0$dose', conf.level = 0.95)
  result_0 <- data.frame(posthoc_0$`data_0$dose`)
  df_tukey_0[1:3,i] <- result_0["p.adj"]
  
  ## time 5
  data_5 <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$time == 5],
                       dose = df_BIT_anova$dose[df_BIT_anova$time == 5])
  lm_5 <- lm(formula = intensity ~ dose, data = data_5)
  df_anova_5[1,i] <- anova(lm_5)$`Pr(>F)`[1]
  a_5 <- aov(data_5$intensity ~ data_5$dose)
  posthoc_5 <- TukeyHSD(x = a_5, 'data_5$dose', conf.level = 0.95)
  result_5 <- data.frame(posthoc_5$`data_5$dose`)
  df_tukey_5[1:3,i] <- result_5["p.adj"]
  
  ## time 24
  data_24 <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$time == 24],
                       dose = df_BIT_anova$dose[df_BIT_anova$time == 24])
  lm_24 <- lm(formula = intensity ~ dose, data = data_24)
  df_anova_24[1,i] <- anova(lm_24)$`Pr(>F)`[1]
  a_24 <- aov(data_24$intensity ~ data_24$dose)
  posthoc_24 <- TukeyHSD(x = a_24, 'data_24$dose', conf.level = 0.95)
  result_24 <- data.frame(posthoc_24$`data_24$dose`)
  df_tukey_24[1:3,i] <- result_24["p.adj"]
}

# Multiple comparison testing using Benjamini Hochberg
# It is necessary to correct for multiple comparisons AFTER running ANOVA
df_anova_0_padj <- as.data.frame(t(df_anova_0))
df_anova_5_padj <- as.data.frame(t(df_anova_5))
df_anova_24_padj <- as.data.frame(t(df_anova_24))
df_anova_0_padj[,1] <- p.adjust(df_anova_0_padj[,1], method = "BH")
df_anova_5_padj[,1] <- p.adjust(df_anova_5_padj[,1], method = "BH")
df_anova_24_padj[,1] <- p.adjust(df_anova_24_padj[,1], method = "BH")

# combine anova and post-hoc results and filter for significant anova results
df_anova_res_0 <- cbind(t(df_tukey_0),df_anova_0_padj)
colnames(df_anova_res_0) <- c("0_0.1-0","0_10-0","0_10-0.1","p.adj")
df_anova_res_0 <- subset(df_anova_res_0, p.adj < 0.05)

df_anova_res_5 <- cbind(t(df_tukey_5),df_anova_5_padj)
colnames(df_anova_res_5) <- c("5_0.1-0","5_10-0","5_10-0.1","p.adj")
df_anova_res_5 <- subset(df_anova_res_5, p.adj < 0.05)

df_anova_res_24 <- cbind(t(df_tukey_24),df_anova_24_padj)
colnames(df_anova_res_24) <- c("24_0.1-0","24_10-0","24_10-0.1","p.adj")
df_anova_res_24 <- subset(df_anova_res_24, p.adj < 0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: merge anova-log2fc #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_FC_0 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_anova_res_0))
df_FC_5 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_anova_res_5))
df_FC_24 <- subset(df_trans_FC_time, row.names(df_trans_FC_time) %in% row.names(df_anova_res_24))

df_BIT_0_0.1v0 <- data.frame(df_FC_0[1], df_anova_res_0[1])
colnames(df_BIT_0_0.1v0) <- c("log2FC","p")
df_BIT_0_0.1v0$condition <- c("0_0.1-0")
df_BIT_0_0.1v0$metab <- rownames(df_BIT_0_0.1v0)

df_BIT_0_10v0 <- data.frame(df_FC_0[2], df_anova_res_0[2])
colnames(df_BIT_0_10v0) <- c("log2FC","p")
df_BIT_0_10v0$condition <- c("0_10-0")
df_BIT_0_10v0$metab <- rownames(df_BIT_0_10v0)

df_BIT_0_10v0.1 <- data.frame(df_FC_0[3], df_anova_res_0[3])
colnames(df_BIT_0_10v0.1) <- c("log2FC","p")
df_BIT_0_10v0.1$condition <- c("0_10-0.1")
df_BIT_0_10v0.1$metab <- rownames(df_BIT_0_10v0.1)

df_BIT_5_0.1v0 <- data.frame(df_FC_5[4], df_anova_res_5[1])
colnames(df_BIT_5_0.1v0) <- c("log2FC","p")
df_BIT_5_0.1v0$condition <- c("5_0.1-0")
df_BIT_5_0.1v0$metab <- rownames(df_BIT_5_0.1v0)

df_BIT_5_10v0 <- data.frame(df_FC_5[5], df_anova_res_5[2])
colnames(df_BIT_5_10v0) <- c("log2FC","p")
df_BIT_5_10v0$condition <- c("5_10-0")
df_BIT_5_10v0$metab <- rownames(df_BIT_5_10v0)

df_BIT_5_10v0.1 <- data.frame(df_FC_5[6], df_anova_res_5[3])
colnames(df_BIT_5_10v0.1) <- c("log2FC","p")
df_BIT_5_10v0.1$condition <- c("5_10-0.1")
df_BIT_5_10v0.1$metab <- rownames(df_BIT_5_10v0.1)

df_BIT_24_0.1v0 <- data.frame(df_FC_24[7], df_anova_res_24[1])
colnames(df_BIT_24_0.1v0) <- c("log2FC","p")
df_BIT_24_0.1v0$condition <- c("24_0.1-0")
df_BIT_24_0.1v0$metab <- rownames(df_BIT_24_0.1v0)

df_BIT_24_10v0 <- data.frame(df_FC_24[8], df_anova_res_24[2])
colnames(df_BIT_24_10v0) <- c("log2FC","p")
df_BIT_24_10v0$condition <- c("24_10-0")
df_BIT_24_10v0$metab <- rownames(df_BIT_24_10v0)

df_BIT_24_10v0.1 <- data.frame(df_FC_24[9], df_anova_res_24[3])
colnames(df_BIT_24_10v0.1) <- c("log2FC","p")
df_BIT_24_10v0.1$condition <- c("24_10-0.1")
df_BIT_24_10v0.1$metab <- rownames(df_BIT_24_10v0.1)

df_BIT_uniStats_0 <- rbind(df_BIT_0_0.1v0, df_BIT_0_10v0, df_BIT_0_10v0.1)
df_BIT_uniStats_0$condition <- factor(df_BIT_uniStats_0$condition,levels=c("0_0.1-0","0_10-0","0_10-0.1"))

df_BIT_uniStats_5 <- rbind(df_BIT_5_0.1v0, df_BIT_5_10v0, df_BIT_5_10v0.1)
df_BIT_uniStats_5$condition <- factor(df_BIT_uniStats_5$condition,levels=c("5_0.1-0","5_10-0","5_10-0.1"))

df_BIT_uniStats_24 <- rbind(df_BIT_24_0.1v0, df_BIT_24_10v0, df_BIT_24_10v0.1)
df_BIT_uniStats_24$condition <- factor(df_BIT_uniStats_24$condition,levels=c("24_0.1-0","24_10-0","24_10-0.1"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: volcano plot for anova #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# inspired from here: https://twbattaglia.github.io/2016/12/17/volcano-plot/

df_BIT_uniStats_0 <- df_BIT_uniStats_0 %>%
  mutate(color = ifelse(df_BIT_uniStats_0$log2FC>2 & df_BIT_uniStats_0$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_0$log2FC<(-2) & df_BIT_uniStats_0$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

df_BIT_uniStats_5 <- df_BIT_uniStats_5 %>%
  mutate(color = ifelse(df_BIT_uniStats_5$log2FC>2 & df_BIT_uniStats_5$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_5$log2FC<(-2) & df_BIT_uniStats_5$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

df_BIT_uniStats_24 <- df_BIT_uniStats_24 %>%
  mutate(color = ifelse(df_BIT_uniStats_24$log2FC>2 & df_BIT_uniStats_24$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_24$log2FC<(-2) & df_BIT_uniStats_24$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

ggplot(df_BIT_uniStats_0, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

ggplot(df_BIT_uniStats_5, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

ggplot(df_BIT_uniStats_24, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### time: heatmap anova #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_0_0.1v0_sig <- filter(df_BIT_0_0.1v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_0_10v0_sig <- filter(df_BIT_0_10v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_0_10v0.1_sig <- filter(df_BIT_0_10v0.1,abs(log2FC) > 2 & p < 0.05)
df_BIT_5_0.1v0_sig <- filter(df_BIT_5_0.1v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_5_10v0_sig <- filter(df_BIT_5_10v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_5_10v0.1_sig <- filter(df_BIT_5_10v0.1,abs(log2FC) > 2 & p < 0.05)
df_BIT_24_0.1v0_sig <- filter(df_BIT_24_0.1v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_24_10v0_sig <- filter(df_BIT_24_10v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_24_10v0.1_sig <- filter(df_BIT_24_10v0.1,abs(log2FC) > 2 & p < 0.05)

sig_metabs_0 <- unique(select(rbind(df_BIT_0_0.1v0_sig, df_BIT_0_10v0_sig, df_BIT_0_10v0.1_sig),metab))
sig_metabs_5 <- unique(select(rbind(df_BIT_5_0.1v0_sig, df_BIT_5_10v0_sig, df_BIT_5_10v0.1_sig),metab))
sig_metabs_24 <- unique(select(rbind(df_BIT_24_0.1v0_sig, df_BIT_24_10v0_sig, df_BIT_24_10v0.1_sig),metab))

heatmap_sig_0 <- sig_metabs_0
heatmap_sig_0$metab<- sort(heatmap_sig_0$metab)
heatmap_sig_5 <- sig_metabs_5
heatmap_sig_5$metab<- sort(heatmap_sig_5$metab)
heatmap_sig_24 <- sig_metabs_24
heatmap_sig_24$metab<- sort(heatmap_sig_24$metab)

heatmap_sig_0$"0_0.1-0" <- rep(0,nrow(heatmap_sig_0))
heatmap_sig_0$"0_10-0" <- rep(0,nrow(heatmap_sig_0))
heatmap_sig_0$"0_10-0.1" <- rep(0,nrow(heatmap_sig_0))
heatmap_sig_5$"5_0.1-0" <- rep(0,nrow(heatmap_sig_5))
heatmap_sig_5$"5_10-0" <- rep(0,nrow(heatmap_sig_5))
heatmap_sig_5$"5_10-0.1" <- rep(0,nrow(heatmap_sig_5))
heatmap_sig_24$"24_0.1-0" <- rep(0,nrow(heatmap_sig_24))
heatmap_sig_24$"24_10-0" <- rep(0,nrow(heatmap_sig_24))
heatmap_sig_24$"24_10-0.1" <- rep(0,nrow(heatmap_sig_24))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
heatmap_sig_0[heatmap_sig_0$metab %in% df_BIT_0_0.1v0$metab,]$"0_0.1-0" = df_BIT_0_0.1v0[df_BIT_0_0.1v0$metab %in% heatmap_sig_0$metab,]$log2FC
heatmap_sig_0[heatmap_sig_0$metab %in% df_BIT_0_10v0$metab,]$"0_10-0" = df_BIT_0_10v0[df_BIT_0_10v0$metab %in% heatmap_sig_0$metab,]$log2FC
heatmap_sig_0[heatmap_sig_0$metab %in% df_BIT_0_10v0.1$metab,]$"0_10-0.1" = df_BIT_0_10v0.1[df_BIT_0_10v0.1$metab %in% heatmap_sig_0$metab,]$log2FC
heatmap_sig_5[heatmap_sig_5$metab %in% df_BIT_5_0.1v0$metab,]$"5_0.1-0" = df_BIT_5_0.1v0[df_BIT_5_0.1v0$metab %in% heatmap_sig_5$metab,]$log2FC
heatmap_sig_5[heatmap_sig_5$metab %in% df_BIT_5_10v0$metab,]$"5_10-0" = df_BIT_5_10v0[df_BIT_5_10v0$metab %in% heatmap_sig_5$metab,]$log2FC
heatmap_sig_5[heatmap_sig_5$metab %in% df_BIT_5_10v0.1$metab,]$"5_10-0.1" = df_BIT_5_10v0.1[df_BIT_5_10v0.1$metab %in% heatmap_sig_5$metab,]$log2FC
heatmap_sig_24[heatmap_sig_24$metab %in% df_BIT_24_0.1v0$metab,]$"24_0.1-0" = df_BIT_24_0.1v0[df_BIT_24_0.1v0$metab %in% heatmap_sig_24$metab,]$log2FC
heatmap_sig_24[heatmap_sig_24$metab %in% df_BIT_24_10v0$metab,]$"24_10-0" = df_BIT_24_10v0[df_BIT_24_10v0$metab %in% heatmap_sig_24$metab,]$log2FC
heatmap_sig_24[heatmap_sig_24$metab %in% df_BIT_24_10v0.1$metab,]$"24_10-0.1" = df_BIT_24_10v0.1[df_BIT_24_10v0.1$metab %in% heatmap_sig_24$metab,]$log2FC

rownames(heatmap_sig_0) <- heatmap_sig_0$metab
heatmap_sig_0 <- as.matrix(select(heatmap_sig_0,-metab))
rownames(heatmap_sig_5) <- heatmap_sig_5$metab
heatmap_sig_5 <- as.matrix(select(heatmap_sig_5,-metab))
rownames(heatmap_sig_24) <- heatmap_sig_24$metab
heatmap_sig_24 <- as.matrix(select(heatmap_sig_24,-metab))

heatmap.2(heatmap_sig_0,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

heatmap.2(heatmap_sig_5,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

heatmap.2(heatmap_sig_24,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: univariate assumptions assessment #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## Normality
# because our sample size is small (8 replicates), statistical normality tests like Shapiro-Wilks are not really appropriate
# we can either assume normality (and run a one-way ANOVA) OR we can perform a non-parametric test (Kruskal-Wallis H Test)

## Sphericity
# use Mauchly's Test for Sphericity, which will test whether variances of the differences between all combinations of related groups are equal

df_mauchly_untreated <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_mauchly_untreated) <- c("p-val")
colnames(df_mauchly_untreated) <- metabs

df_mauchly_persister <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_mauchly_persister) <- c("p-val")
colnames(df_mauchly_persister) <- metabs

df_mauchly_dead <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_mauchly_dead) <- c("pval")
colnames(df_mauchly_dead) <- metabs

for (i in 1:length(metabs)) {
  
  timeLevels <- c(1,2,3)
  timeFactor <- as.factor(timeLevels)
  timeFrame <- data.frame(timeFactor)
  
  ## untreated
  untreatedBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 0],
                    df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 5],
                    df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 24])
  untreatedModel <- lm(untreatedBind ~ 1)
  untreatedAnalysis <- Anova(untreatedModel, idata = timeFrame, idesign = ~timeFactor)
  df_mauchly_untreated[i] <- unlist(summary(untreatedAnalysis)[6])[2]
  
  ## persister
  persisterBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 0],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 5],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 24])
  persisterModel <- lm(persisterBind ~ 1)
  persisterAnalysis <- Anova(persisterModel, idata = timeFrame, idesign = ~timeFactor)
  df_mauchly_persister[i] <- unlist(summary(persisterAnalysis)[6])[2]
  
  ## dead
  deadBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 0],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 5],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 24])
  deadModel <- lm(deadBind ~ 1)
  deadAnalysis <- Anova(deadModel, idata = timeFrame, idesign = ~timeFactor)
  df_mauchly_dead[i] <- unlist(summary(deadAnalysis)[6])[2]
}

length(which(df_mauchly_untreated < 0.05)) # 106
length(which(df_mauchly_persister < 0.05)) # 83
length(which(df_mauchly_dead < 0.05)) # 74

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: repeated measures one-way ANOVA with Greenhouse-Geiser correction and Games-Howell post-hoc test #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reference: Repeated Measures ANOVA: Introduction to Statistics Using R (Pyschology 9041B), Paul Gribble
# https://www.gribblelab.org/stats/notes/RepeatedMeasuresANOVA.pdf

df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## One-way ANOVA for the different time-points
#df_BIT_anova <- unite(df_BIT_uni, group, c("dose"))
df_BIT_greenhouse <- df_BIT_uni
df_BIT_greenhouse <- select(df_BIT_greenhouse, -techRep)

# create a dataframe to save the sample ANOVA results
df_greenhouse_untreated <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_greenhouse_untreated) <- c("untreated")
colnames(df_greenhouse_untreated) <- metabs

df_greenhouse_persister <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_greenhouse_persister) <- c("persister")
colnames(df_greenhouse_persister) <- metabs

df_greenhouse_dead <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_greenhouse_dead) <- c("dead")
colnames(df_greenhouse_dead) <- metabs

df_games_untreated <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_untreated) <- c("untreated_24-0","untreated_5-0","untreated_5-24")
colnames(df_games_untreated) <- metabs

df_games_persister <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_persister) <- c("persister_24-0","persister_5-0","persister_5-24")
colnames(df_games_persister) <- metabs

df_games_dead <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_games_dead) <- c("dead_24-0","dead_5-0","dead_5-24")
colnames(df_games_dead) <- metabs

# run ANOVA for all the untreated, persister, and dead conditions
for (i in 1:length(metabs)) {
  
  timeLevels <- c(1,2,3)
  timeFactor <- as.factor(timeLevels)
  timeFrame <- data.frame(timeFactor)
  
  ## untreated
  untreatedBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 0],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 5],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0 & df_BIT_uni$time == 24])
  untreatedModel <- lm(untreatedBind ~ 1)
  untreatedAnalysis <- Anova(untreatedModel, idata = timeFrame, idesign = ~timeFactor)
  df_greenhouse_untreated[i] <- unlist(summary(untreatedAnalysis)[5])[2]
  # games-howell post-hoc
  data_untreated <- data.frame(intensity = df_BIT_greenhouse[,metabs[i]][df_BIT_greenhouse$dose == 0],
                               time = df_BIT_greenhouse$time[df_BIT_greenhouse$dose == 0])
  posthoc_untreated <- posthocTGH(data_untreated$intensity, data_untreated$time, method=c("games-howell"), digits = 5)
  df_games_untreated[1,i] <- unlist(posthoc_untreated$output[2])[16]
  df_games_untreated[2,i] <- unlist(posthoc_untreated$output[2])[17]
  df_games_untreated[3,i] <- unlist(posthoc_untreated$output[2])[18]
  
  ## persister
  persisterBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 0],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 5],
                         df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==0.1 & df_BIT_uni$time == 24])
  persisterModel <- lm(persisterBind ~ 1)
  persisterAnalysis <- Anova(persisterModel, idata = timeFrame, idesign = ~timeFactor)
  df_greenhouse_persister[i] <- unlist(summary(persisterAnalysis)[5])[2]
  # games-howell post-hoc
  data_persister <- data.frame(intensity = df_BIT_greenhouse[,metabs[i]][df_BIT_greenhouse$dose == 0.1],
                               time = df_BIT_greenhouse$time[df_BIT_greenhouse$dose == 0.1])
  posthoc_persister <- posthocTGH(data_persister$intensity, data_persister$time, method=c("games-howell"), digits = 5)
  df_games_persister[1,i] <- unlist(posthoc_persister$output[2])[16]
  df_games_persister[2,i] <- unlist(posthoc_persister$output[2])[17]
  df_games_persister[3,i] <- unlist(posthoc_persister$output[2])[18]
  
  ## dead
  deadBind <- cbind(df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 0],
                    df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 5],
                    df_BIT_uni[,metabs[i]][df_BIT_uni$dose ==10 & df_BIT_uni$time == 24])
  deadModel <- lm(deadBind ~ 1)
  deadAnalysis <- Anova(deadModel, idata = timeFrame, idesign = ~timeFactor)
  df_greenhouse_dead[i] <- unlist(summary(deadAnalysis)[5])[2]
  # games-howell post-hoc
  data_dead <- data.frame(intensity = df_BIT_greenhouse[,metabs[i]][df_BIT_greenhouse$dose == 10],
                               time = df_BIT_greenhouse$time[df_BIT_greenhouse$dose == 10])
  posthoc_dead <- posthocTGH(data_dead$intensity, data_dead$time, method=c("games-howell"), digits = 5)
  df_games_dead[1,i] <- unlist(posthoc_dead$output[2])[16]
  df_games_dead[2,i] <- unlist(posthoc_dead$output[2])[17]
  df_games_dead[3,i] <- unlist(posthoc_dead$output[2])[18]
}

# Multiple comparison testing using Benjamini Hochberg
# It is necessary to correct for multiple comparisons AFTER running ANOVA
df_greenhouse_untreated_padj <- as.data.frame(t(df_greenhouse_untreated))
df_greenhouse_persister_padj <- as.data.frame(t(df_greenhouse_persister))
df_greenhouse_dead_padj <- as.data.frame(t(df_greenhouse_dead))
df_greenhouse_untreated_padj[,1] <- p.adjust(df_greenhouse_untreated_padj[,1], method = "BH")
df_greenhouse_persister_padj[,1] <- p.adjust(df_greenhouse_persister_padj[,1], method = "BH")
df_greenhouse_dead_padj[,1] <- p.adjust(df_greenhouse_dead_padj[,1], method = "BH")

# combine anova and post-hoc results and filter for significant anova results
df_greenhouse_res_untreated <- cbind(t(df_games_untreated),df_greenhouse_untreated_padj)
colnames(df_greenhouse_res_untreated) <- c("untreated_24-0","untreated_5-0","untreated_5-24","p.adj")
df_greenhouse_res_untreated <- subset(df_greenhouse_res_untreated, p.adj < 0.05)

df_greenhouse_res_persister <- cbind(t(df_games_persister),df_greenhouse_persister_padj)
colnames(df_greenhouse_res_persister) <- c("persister_24-0","persister_5-0","persister_5-24","p.adj")
df_greenhouse_res_persister <- subset(df_greenhouse_res_persister, p.adj < 0.05)

df_greenhouse_res_dead <- cbind(t(df_games_dead),df_greenhouse_dead_padj)
colnames(df_greenhouse_res_dead) <- c("dead_24-0","dead_5-0","dead_5-24","p.adj")
df_greenhouse_res_dead <- subset(df_greenhouse_res_dead, p.adj < 0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: merge anova-log2fc #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_FC_untreated <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_greenhouse_res_untreated))
df_FC_persister <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_greenhouse_res_persister))
df_FC_dead <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_greenhouse_res_dead))

df_BIT_untreated_24v0 <- data.frame(df_FC_untreated[1], df_greenhouse_res_untreated[1])
colnames(df_BIT_untreated_24v0) <- c("log2FC","p")
df_BIT_untreated_24v0$condition <- c("untreated_24-0")
df_BIT_untreated_24v0$metab <- rownames(df_BIT_untreated_24v0)

df_BIT_untreated_5v0 <- data.frame(df_FC_untreated[2], df_greenhouse_res_untreated[2])
colnames(df_BIT_untreated_5v0) <- c("log2FC","p")
df_BIT_untreated_5v0$condition <- c("untreated_5-0")
df_BIT_untreated_5v0$metab <- rownames(df_BIT_untreated_5v0)

df_BIT_untreated_5v24 <- data.frame(df_FC_untreated[3], df_greenhouse_res_untreated[3])
colnames(df_BIT_untreated_5v24) <- c("log2FC","p")
df_BIT_untreated_5v24$condition <- c("untreated_5-24")
df_BIT_untreated_5v24$metab <- rownames(df_BIT_untreated_5v24)

df_BIT_persister_24v0 <- data.frame(df_FC_persister[4], df_greenhouse_res_persister[1])
colnames(df_BIT_persister_24v0) <- c("log2FC","p")
df_BIT_persister_24v0$condition <- c("persister_24-0")
df_BIT_persister_24v0$metab <- rownames(df_BIT_persister_24v0)

df_BIT_persister_5v0 <- data.frame(df_FC_persister[5], df_greenhouse_res_persister[2])
colnames(df_BIT_persister_5v0) <- c("log2FC","p")
df_BIT_persister_5v0$condition <- c("persister_5-0")
df_BIT_persister_5v0$metab <- rownames(df_BIT_persister_5v0)

df_BIT_persister_5v24 <- data.frame(df_FC_persister[6], df_greenhouse_res_persister[3])
colnames(df_BIT_persister_5v24) <- c("log2FC","p")
df_BIT_persister_5v24$condition <- c("persister_5-24")
df_BIT_persister_5v24$metab <- rownames(df_BIT_persister_5v24)

df_BIT_dead_24v0 <- data.frame(df_FC_dead[7], df_greenhouse_res_dead[1])
colnames(df_BIT_dead_24v0) <- c("log2FC","p")
df_BIT_dead_24v0$condition <- c("dead_24-0")
df_BIT_dead_24v0$metab <- rownames(df_BIT_dead_24v0)

df_BIT_dead_5v0 <- data.frame(df_FC_dead[8], df_greenhouse_res_dead[2])
colnames(df_BIT_dead_5v0) <- c("log2FC","p")
df_BIT_dead_5v0$condition <- c("dead_5-0")
df_BIT_dead_5v0$metab <- rownames(df_BIT_dead_5v0)

df_BIT_dead_5v24 <- data.frame(df_FC_dead[9], df_greenhouse_res_dead[3])
colnames(df_BIT_dead_5v24) <- c("log2FC","p")
df_BIT_dead_5v24$condition <- c("dead_5-24")
df_BIT_dead_5v24$metab <- rownames(df_BIT_dead_5v24)

df_BIT_uniStats_untreated <- rbind(df_BIT_untreated_24v0, df_BIT_untreated_5v0, df_BIT_untreated_5v24)
df_BIT_uniStats_untreated$condition <- factor(df_BIT_uniStats_untreated$condition,levels=c("untreated_24-0","untreated_5-0","untreated_5-24"))

df_BIT_uniStats_persister <- rbind(df_BIT_persister_24v0, df_BIT_persister_5v0, df_BIT_persister_5v24)
df_BIT_uniStats_persister$condition <- factor(df_BIT_uniStats_persister$condition,levels=c("persister_24-0","persister_5-0","persister_5-24"))

df_BIT_uniStats_dead <- rbind(df_BIT_dead_24v0, df_BIT_dead_5v0, df_BIT_dead_5v24)
df_BIT_uniStats_dead$condition <- factor(df_BIT_uniStats_dead$condition,levels=c("dead_24-0","dead_5-0","dead_5-24"))

write.csv(df_BIT_uniStats_untreated, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_untreated.csv'))
write.csv(df_BIT_uniStats_persister, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_persister.csv'))
write.csv(df_BIT_uniStats_dead, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_dead.csv'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: repeated measures one-way ANOVA with Tukey post-hoc test #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reference: Repeated Measures ANOVA: Introduction to Statistics Using R (Pyschology 9041B), Paul Gribble
# https://www.gribblelab.org/stats/notes/RepeatedMeasuresANOVA.pdf

df_BIT_uni <- as.data.frame(df_BIT_trans)

# create a new variable, sample, that is the ID names
df_BIT_uni$sample <- rownames(df_BIT_uni)
# remove row names
rownames(df_BIT_uni) <- c()
# separate the IDs for ease of calling
df_BIT_uni <- separate(df_BIT_uni, sample, into=c("dose","time","bioRep","techRep"), sep = "_")

# obtain metabolite IDs
metabs <- colnames(df_BIT_uni)[1:336]

## One-way ANOVA for the different time-points
#df_BIT_anova <- unite(df_BIT_uni, group, c("dose"))
df_BIT_anova <- df_BIT_uni
df_BIT_anova <- select(df_BIT_anova, -techRep)

# create a dataframe to save the sample ANOVA results
df_anova_untreated <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_untreated) <- c("untreated")
colnames(df_anova_untreated) <- metabs

df_anova_persister <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_persister) <- c("persister")
colnames(df_anova_persister) <- metabs

df_anova_dead <- data.frame(matrix(NA, nrow = 1, ncol = length(metabs)))
row.names(df_anova_dead) <- c("dead")
colnames(df_anova_dead) <- metabs

df_tukey_untreated <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_untreated) <- c("untreated_24-0","untreated_5-0","untreated_5-24")
colnames(df_tukey_untreated) <- metabs

df_tukey_persister <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_persister) <- c("persister_24-0","persister_5-0","persister_5-24")
colnames(df_tukey_persister) <- metabs

df_tukey_dead <- data.frame(matrix(NA, nrow = 3, ncol = length(metabs)))
row.names(df_tukey_dead) <- c("dead_24-0","dead_5-0","dead_5-24")
colnames(df_tukey_dead) <- metabs

# run ANOVA for all the untreated, persister, and dead conditions
for (i in 1:length(metabs)) {
  
  ## untreated
  data_untreated <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$dose == 0],
                               time = df_BIT_anova[,338][df_BIT_anova$dose == 0],
                               bioRep = df_BIT_anova[,339][df_BIT_anova$dose == 0])
  data_untreated$time <- as.factor(data_untreated$time)
  data_untreated_anova <- aov(intensity ~ time + Error(bioRep/time), data = data_untreated)
  df_anova_untreated[1,i]<- unlist(summary(data_untreated_anova)[2])[9]
  data_untreated_lme <- lme(intensity ~ time, random = ~1|bioRep/time, data = data_untreated)
  #anova(data_untreated_lme)
  df_tukey_untreated[1,i] <- unlist(summary(glht(data_untreated_lme, linfct = mcp(time = "Tukey")))[10])[12]
  df_tukey_untreated[2,i] <- unlist(summary(glht(data_untreated_lme, linfct = mcp(time = "Tukey")))[10])[13]
  df_tukey_untreated[3,i] <- unlist(summary(glht(data_untreated_lme, linfct = mcp(time = "Tukey")))[10])[14]
  
  ## persister
  data_persister <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$dose == 0.1],
                               time = df_BIT_anova[,338][df_BIT_anova$dose == 0.1],
                               bioRep = df_BIT_anova[,339][df_BIT_anova$dose == 0.1])
  data_persister$time <- as.factor(data_persister$time)
  data_persister_anova <- aov(intensity ~ time + Error(bioRep/time), data = data_persister)
  df_anova_persister[1,i]<- unlist(summary(data_persister_anova)[2])[9]
  data_persister_lme <- lme(intensity ~ time, random = ~1|bioRep/time, data = data_persister)
  #anova(data_persister_lme)
  df_tukey_persister[1,i] <- unlist(summary(glht(data_persister_lme, linfct = mcp(time = "Tukey")))[10])[12]
  df_tukey_persister[2,i] <- unlist(summary(glht(data_persister_lme, linfct = mcp(time = "Tukey")))[10])[13]
  df_tukey_persister[3,i] <- unlist(summary(glht(data_persister_lme, linfct = mcp(time = "Tukey")))[10])[14]
  
  ## dead
  data_dead <- data.frame(intensity = df_BIT_anova[,metabs[i]][df_BIT_anova$dose == 10],
                               time = df_BIT_anova[,338][df_BIT_anova$dose == 10],
                               bioRep = df_BIT_anova[,339][df_BIT_anova$dose == 10])
  data_dead$time <- as.factor(data_dead$time)
  data_dead_anova <- aov(intensity ~ time + Error(bioRep/time), data = data_dead)
  df_anova_dead[1,i]<- unlist(summary(data_dead_anova)[2])[9]
  data_dead_lme <- lme(intensity ~ time, random = ~1|bioRep/time, data = data_dead)
  #anova(data_dead_lme)
  df_tukey_dead[1,i] <- unlist(summary(glht(data_dead_lme, linfct = mcp(time = "Tukey")))[10])[12]
  df_tukey_dead[2,i] <- unlist(summary(glht(data_dead_lme, linfct = mcp(time = "Tukey")))[10])[13]
  df_tukey_dead[3,i] <- unlist(summary(glht(data_dead_lme, linfct = mcp(time = "Tukey")))[10])[14]

}

# Multiple comparison testing using Benjamini Hochberg
# It is necessary to correct for multiple comparisons AFTER running ANOVA
df_anova_untreated_padj <- as.data.frame(t(df_anova_untreated))
df_anova_persister_padj <- as.data.frame(t(df_anova_persister))
df_anova_dead_padj <- as.data.frame(t(df_anova_dead))
df_anova_untreated_padj[,1] <- p.adjust(df_anova_untreated_padj[,1], method = "BH")
df_anova_persister_padj[,1] <- p.adjust(df_anova_persister_padj[,1], method = "BH")
df_anova_dead_padj[,1] <- p.adjust(df_anova_dead_padj[,1], method = "BH")

# combine anova and post-hoc results and filter for significant anova results
df_anova_res_untreated <- cbind(t(df_tukey_untreated),df_anova_untreated_padj)
colnames(df_anova_res_untreated) <- c("untreated_24-0","untreated_5-0","untreated_5-24","p.adj")
df_anova_res_untreated <- subset(df_anova_res_untreated, p.adj < 0.05)

df_anova_res_persister <- cbind(t(df_tukey_persister),df_anova_persister_padj)
colnames(df_anova_res_persister) <- c("persister_24-0","persister_5-0","persister_5-24","p.adj")
df_anova_res_persister <- subset(df_anova_res_persister, p.adj < 0.05)

df_anova_res_dead <- cbind(t(df_tukey_dead),df_anova_dead_padj)
colnames(df_anova_res_dead) <- c("dead_24-0","dead_5-0","dead_5-24","p.adj")
df_anova_res_dead <- subset(df_anova_res_dead, p.adj < 0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: merge anova-log2fc #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_FC_untreated <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_anova_res_untreated))
df_FC_persister <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_anova_res_persister))
df_FC_dead <- subset(df_trans_FC_sample, row.names(df_trans_FC_sample) %in% row.names(df_anova_res_dead))

df_BIT_untreated_24v0 <- data.frame(df_FC_untreated[1], df_anova_res_untreated[1])
colnames(df_BIT_untreated_24v0) <- c("log2FC","p")
df_BIT_untreated_24v0$condition <- c("untreated_24-0")
df_BIT_untreated_24v0$metab <- rownames(df_BIT_untreated_24v0)

df_BIT_untreated_5v0 <- data.frame(df_FC_untreated[2], df_anova_res_untreated[2])
colnames(df_BIT_untreated_5v0) <- c("log2FC","p")
df_BIT_untreated_5v0$condition <- c("untreated_5-0")
df_BIT_untreated_5v0$metab <- rownames(df_BIT_untreated_5v0)

df_BIT_untreated_5v24 <- data.frame(df_FC_untreated[3], df_anova_res_untreated[3])
colnames(df_BIT_untreated_5v24) <- c("log2FC","p")
df_BIT_untreated_5v24$condition <- c("untreated_5-24")
df_BIT_untreated_5v24$metab <- rownames(df_BIT_untreated_5v24)

df_BIT_persister_24v0 <- data.frame(df_FC_persister[4], df_anova_res_persister[1])
colnames(df_BIT_persister_24v0) <- c("log2FC","p")
df_BIT_persister_24v0$condition <- c("persister_24-0")
df_BIT_persister_24v0$metab <- rownames(df_BIT_persister_24v0)

df_BIT_persister_5v0 <- data.frame(df_FC_persister[5], df_anova_res_persister[2])
colnames(df_BIT_persister_5v0) <- c("log2FC","p")
df_BIT_persister_5v0$condition <- c("persister_5-0")
df_BIT_persister_5v0$metab <- rownames(df_BIT_persister_5v0)

df_BIT_persister_5v24 <- data.frame(df_FC_persister[6], df_anova_res_persister[3])
colnames(df_BIT_persister_5v24) <- c("log2FC","p")
df_BIT_persister_5v24$condition <- c("persister_5-24")
df_BIT_persister_5v24$metab <- rownames(df_BIT_persister_5v24)

df_BIT_dead_24v0 <- data.frame(df_FC_dead[7], df_anova_res_dead[1])
colnames(df_BIT_dead_24v0) <- c("log2FC","p")
df_BIT_dead_24v0$condition <- c("dead_24-0")
df_BIT_dead_24v0$metab <- rownames(df_BIT_dead_24v0)

df_BIT_dead_5v0 <- data.frame(df_FC_dead[8], df_anova_res_dead[2])
colnames(df_BIT_dead_5v0) <- c("log2FC","p")
df_BIT_dead_5v0$condition <- c("dead_5-0")
df_BIT_dead_5v0$metab <- rownames(df_BIT_dead_5v0)

df_BIT_dead_5v24 <- data.frame(df_FC_dead[9], df_anova_res_dead[3])
colnames(df_BIT_dead_5v24) <- c("log2FC","p")
df_BIT_dead_5v24$condition <- c("dead_5-24")
df_BIT_dead_5v24$metab <- rownames(df_BIT_dead_5v24)

df_BIT_uniStats_untreated <- rbind(df_BIT_untreated_24v0, df_BIT_untreated_5v0, df_BIT_untreated_5v24)
df_BIT_uniStats_untreated$condition <- factor(df_BIT_uniStats_untreated$condition,levels=c("untreated_24-0","untreated_5-0","untreated_5-24"))

df_BIT_uniStats_persister <- rbind(df_BIT_persister_24v0, df_BIT_persister_5v0, df_BIT_persister_5v24)
df_BIT_uniStats_persister$condition <- factor(df_BIT_uniStats_persister$condition,levels=c("persister_24-0","persister_5-0","persister_5-24"))

df_BIT_uniStats_dead <- rbind(df_BIT_dead_24v0, df_BIT_dead_5v0, df_BIT_dead_5v24)
df_BIT_uniStats_dead$condition <- factor(df_BIT_uniStats_dead$condition,levels=c("dead_24-0","dead_5-0","dead_5-24"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: volcano plot for anova #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# inspired from here: https://twbattaglia.github.io/2016/12/17/volcano-plot/

df_BIT_uniStats_untreated <- df_BIT_uniStats_untreated %>%
  mutate(color = ifelse(df_BIT_uniStats_untreated$log2FC>2 & df_BIT_uniStats_untreated$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_untreated$log2FC<(-2) & df_BIT_uniStats_untreated$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

df_BIT_uniStats_persister <- df_BIT_uniStats_persister %>%
  mutate(color = ifelse(df_BIT_uniStats_persister$log2FC>2 & df_BIT_uniStats_persister$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_persister$log2FC<(-2) & df_BIT_uniStats_persister$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

df_BIT_uniStats_dead <- df_BIT_uniStats_dead %>%
  mutate(color = ifelse(df_BIT_uniStats_dead$log2FC>2 & df_BIT_uniStats_dead$p < 0.05,
                        yes = "num",
                        no = ifelse(df_BIT_uniStats_dead$log2FC<(-2) & df_BIT_uniStats_dead$p < 0.05,
                                    yes = "denom",
                                    no = "none")))

ggplot(df_BIT_uniStats_untreated, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

ggplot(df_BIT_uniStats_persister, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

ggplot(df_BIT_uniStats_dead, aes(x=log2FC, y=-log10(p))) +
  geom_point(aes(color = factor(color)),alpha=0.4, size=1.75, na.rm=T) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"),
        strip.text=element_text(face = "bold",size=14),
        strip.background = element_rect(fill="white")) +
  xlab("log2(Fold Change)") +
  scale_color_manual(values = c("num" = myPalette[10],
                                "denom" = myPalette[9],
                                "none" = myPalette[2]),guide=F) +
  #geom_vline(xintercept = 0, colour = "black") + # add line at 0
  #geom_hline(yintercept = 1.3, colour = "black") +
  facet_wrap(~condition,nrow=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### sample: heatmap anova #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_BIT_untreated_24v0_sig <- filter(df_BIT_untreated_24v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_untreated_5v0_sig <- filter(df_BIT_untreated_5v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_untreated_5v24_sig <- filter(df_BIT_untreated_5v24,abs(log2FC) > 2 & p < 0.05)
df_BIT_persister_24v0_sig <- filter(df_BIT_persister_24v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_persister_5v0_sig <- filter(df_BIT_persister_5v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_persister_5v24_sig <- filter(df_BIT_persister_5v24,abs(log2FC) > 2 & p < 0.05)
df_BIT_dead_24v0_sig <- filter(df_BIT_dead_24v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_dead_5v0_sig <- filter(df_BIT_dead_5v0,abs(log2FC) > 2 & p < 0.05)
df_BIT_dead_5v24_sig <- filter(df_BIT_dead_5v24,abs(log2FC) > 2 & p < 0.05)

sig_metabs_untreated <- unique(dplyr::select(rbind(df_BIT_untreated_24v0_sig, df_BIT_untreated_5v0_sig, df_BIT_untreated_5v24_sig),metab))
sig_metabs_persister <- unique(dplyr::select(rbind(df_BIT_persister_24v0_sig, df_BIT_persister_5v0_sig, df_BIT_persister_5v24_sig),metab))
sig_metabs_dead <- unique(dplyr::select(rbind(df_BIT_dead_24v0_sig, df_BIT_dead_5v0_sig, df_BIT_dead_5v24_sig),metab))

heatmap_sig_untreated <- sig_metabs_untreated
heatmap_sig_untreated$metab<- sort(heatmap_sig_untreated$metab)
heatmap_sig_persister <- sig_metabs_persister
heatmap_sig_persister$metab<- sort(heatmap_sig_persister$metab)
heatmap_sig_dead <- sig_metabs_dead
heatmap_sig_dead$metab<- sort(heatmap_sig_dead$metab)

heatmap_sig_untreated$"untreated_24-0" <- rep(0,nrow(heatmap_sig_untreated))
heatmap_sig_untreated$"untreated_5-0" <- rep(0,nrow(heatmap_sig_untreated))
heatmap_sig_untreated$"untreated_5-24" <- rep(0,nrow(heatmap_sig_untreated))
heatmap_sig_persister$"persister_24-0" <- rep(0,nrow(heatmap_sig_persister))
heatmap_sig_persister$"persister_5-0" <- rep(0,nrow(heatmap_sig_persister))
heatmap_sig_persister$"persister_5-24" <- rep(0,nrow(heatmap_sig_persister))
heatmap_sig_dead$"dead_24-0" <- rep(0,nrow(heatmap_sig_dead))
heatmap_sig_dead$"dead_5-0" <- rep(0,nrow(heatmap_sig_dead))
heatmap_sig_dead$"dead_5-24" <- rep(0,nrow(heatmap_sig_dead))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
heatmap_sig_untreated[heatmap_sig_untreated$metab %in% df_BIT_untreated_24v0$metab,]$"untreated_24-0" = df_BIT_untreated_24v0[df_BIT_untreated_24v0$metab %in% heatmap_sig_untreated$metab,]$log2FC
heatmap_sig_untreated[heatmap_sig_untreated$metab %in% df_BIT_untreated_5v0$metab,]$"untreated_5-0" = df_BIT_untreated_5v0[df_BIT_untreated_5v0$metab %in% heatmap_sig_untreated$metab,]$log2FC
heatmap_sig_untreated[heatmap_sig_untreated$metab %in% df_BIT_untreated_5v24$metab,]$"untreated_5-24" = df_BIT_untreated_5v24[df_BIT_untreated_5v24$metab %in% heatmap_sig_untreated$metab,]$log2FC
heatmap_sig_persister[heatmap_sig_persister$metab %in% df_BIT_persister_24v0$metab,]$"persister_24-0" = df_BIT_persister_24v0[df_BIT_persister_24v0$metab %in% heatmap_sig_persister$metab,]$log2FC
heatmap_sig_persister[heatmap_sig_persister$metab %in% df_BIT_persister_5v0$metab,]$"persister_5-0" = df_BIT_persister_5v0[df_BIT_persister_5v0$metab %in% heatmap_sig_persister$metab,]$log2FC
heatmap_sig_persister[heatmap_sig_persister$metab %in% df_BIT_persister_5v24$metab,]$"persister_5-24" = df_BIT_persister_5v24[df_BIT_persister_5v24$metab %in% heatmap_sig_persister$metab,]$log2FC
heatmap_sig_dead[heatmap_sig_dead$metab %in% df_BIT_dead_24v0$metab,]$"dead_24-0" = df_BIT_dead_24v0[df_BIT_dead_24v0$metab %in% heatmap_sig_dead$metab,]$log2FC
heatmap_sig_dead[heatmap_sig_dead$metab %in% df_BIT_dead_5v0$metab,]$"dead_5-0" = df_BIT_dead_5v0[df_BIT_dead_5v0$metab %in% heatmap_sig_dead$metab,]$log2FC
heatmap_sig_dead[heatmap_sig_dead$metab %in% df_BIT_dead_5v24$metab,]$"dead_5-24" = df_BIT_dead_5v24[df_BIT_dead_5v24$metab %in% heatmap_sig_dead$metab,]$log2FC

rownames(heatmap_sig_untreated) <- heatmap_sig_untreated$metab
heatmap_sig_untreated <- as.matrix(dplyr::select(heatmap_sig_untreated,-metab))
rownames(heatmap_sig_persister) <- heatmap_sig_persister$metab
heatmap_sig_persister <- as.matrix(dplyr::select(heatmap_sig_persister,-metab))
rownames(heatmap_sig_dead) <- heatmap_sig_dead$metab
heatmap_sig_dead <- as.matrix(dplyr::select(heatmap_sig_dead,-metab))

heatmap.2(heatmap_sig_untreated,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

heatmap.2(heatmap_sig_persister,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

heatmap.2(heatmap_sig_dead,
          col = colorRampPalette(c(myPalette[9],"white",myPalette[10])) (n = 171),
          #scale = "none",
          #ColSideColors = c(myPalette[3],myPalette[11],myPalette[10],myPalette[10],myPalette[4],myPalette[11]),
          key = TRUE,
          key.title = NA,
          key.xlab = c("log2FC"),
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexCol = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### anova data output #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE: running this section requires that you have run the heatmap sections

output_sig_0 <- sig_metabs_0
output_sig_0$metab<- sort(output_sig_0$metab)
output_sig_5 <- sig_metabs_5
output_sig_5$metab<- sort(output_sig_5$metab)
output_sig_24 <- sig_metabs_24
output_sig_24$metab<- sort(output_sig_24$metab)
output_sig_untreated <- sig_metabs_untreated
output_sig_untreated$metab<- sort(output_sig_untreated$metab)
output_sig_persister <- sig_metabs_persister
output_sig_persister$metab<- sort(output_sig_persister$metab)
output_sig_dead <- sig_metabs_dead
output_sig_dead$metab<- sort(output_sig_dead$metab)

output_sig_0$"0_0.1-0_log2FC" <- rep(0,nrow(output_sig_0))
output_sig_0$"0_0.1-0_p" <- rep(0,nrow(output_sig_0))
output_sig_0$"0_10-0_log2FC" <- rep(0,nrow(output_sig_0))
output_sig_0$"0_10-0_p" <- rep(0,nrow(output_sig_0))
output_sig_0$"0_10-0.1_log2FC" <- rep(0,nrow(output_sig_0))
output_sig_0$"0_10-0.1_p" <- rep(0,nrow(output_sig_0))

output_sig_5$"5_0.1-0_log2FC" <- rep(0,nrow(output_sig_5))
output_sig_5$"5_0.1-0_p" <- rep(0,nrow(output_sig_5))
output_sig_5$"5_10-0_log2FC" <- rep(0,nrow(output_sig_5))
output_sig_5$"5_10-0_p" <- rep(0,nrow(output_sig_5))
output_sig_5$"5_10-0.1_log2FC" <- rep(0,nrow(output_sig_5))
output_sig_5$"5_10-0.1_p" <- rep(0,nrow(output_sig_5))

output_sig_24$"24_0.1-0_log2FC" <- rep(0,nrow(output_sig_24))
output_sig_24$"24_0.1-0_p" <- rep(0,nrow(output_sig_24))
output_sig_24$"24_10-0_log2FC" <- rep(0,nrow(output_sig_24))
output_sig_24$"24_10-0_p" <- rep(0,nrow(output_sig_24))
output_sig_24$"24_10-0.1_log2FC" <- rep(0,nrow(output_sig_24))
output_sig_24$"24_10-0.1_p" <- rep(0,nrow(output_sig_24))

output_sig_untreated$"untreated_24-0_log2FC" <- rep(0,nrow(output_sig_untreated))
output_sig_untreated$"untreated_24-0_p" <- rep(0,nrow(output_sig_untreated))
output_sig_untreated$"untreated_5-0_log2FC" <- rep(0,nrow(output_sig_untreated))
output_sig_untreated$"untreated_5-0_p" <- rep(0,nrow(output_sig_untreated))
output_sig_untreated$"untreated_5-24_log2FC" <- rep(0,nrow(output_sig_untreated))
output_sig_untreated$"untreated_5-24_p" <- rep(0,nrow(output_sig_untreated))

output_sig_persister$"persister_24-0_log2FC" <- rep(0,nrow(output_sig_persister))
output_sig_persister$"persister_24-0_p" <- rep(0,nrow(output_sig_persister))
output_sig_persister$"persister_5-0_log2FC" <- rep(0,nrow(output_sig_persister))
output_sig_persister$"persister_5-0_p" <- rep(0,nrow(output_sig_persister))
output_sig_persister$"persister_5-24_log2FC" <- rep(0,nrow(output_sig_persister))
output_sig_persister$"persister_5-24_p" <- rep(0,nrow(output_sig_persister))

output_sig_dead$"dead_24-0_log2FC" <- rep(0,nrow(output_sig_dead))
output_sig_dead$"dead_24-0_p" <- rep(0,nrow(output_sig_dead))
output_sig_dead$"dead_5-0_log2FC" <- rep(0,nrow(output_sig_dead))
output_sig_dead$"dead_5-0_p" <- rep(0,nrow(output_sig_dead))
output_sig_dead$"dead_5-24_log2FC" <- rep(0,nrow(output_sig_dead))
output_sig_dead$"dead_5-24_p" <- rep(0,nrow(output_sig_dead))

# see piping help here: https://stackoverflow.com/questions/13774773/check-whether-value-exist-in-one-data-frame-or-not
output_sig_0[output_sig_0$metab %in% df_BIT_0_0.1v0$metab,]$"0_0.1-0_log2FC" = df_BIT_0_0.1v0[df_BIT_0_0.1v0$metab %in% output_sig_0$metab,]$log2FC
output_sig_0[output_sig_0$metab %in% df_BIT_0_10v0$metab,]$"0_10-0_log2FC" = df_BIT_0_10v0[df_BIT_0_10v0$metab %in% output_sig_0$metab,]$log2FC
output_sig_0[output_sig_0$metab %in% df_BIT_0_10v0.1$metab,]$"0_10-0.1_log2FC" = df_BIT_0_10v0.1[df_BIT_0_10v0.1$metab %in% output_sig_0$metab,]$log2FC
output_sig_5[output_sig_5$metab %in% df_BIT_5_0.1v0$metab,]$"5_0.1-0_log2FC" = df_BIT_5_0.1v0[df_BIT_5_0.1v0$metab %in% output_sig_5$metab,]$log2FC
output_sig_5[output_sig_5$metab %in% df_BIT_5_10v0$metab,]$"5_10-0_log2FC" = df_BIT_5_10v0[df_BIT_5_10v0$metab %in% output_sig_5$metab,]$log2FC
output_sig_5[output_sig_5$metab %in% df_BIT_5_10v0.1$metab,]$"5_10-0.1_log2FC" = df_BIT_5_10v0.1[df_BIT_5_10v0.1$metab %in% output_sig_5$metab,]$log2FC
output_sig_24[output_sig_24$metab %in% df_BIT_24_0.1v0$metab,]$"24_0.1-0_log2FC" = df_BIT_24_0.1v0[df_BIT_24_0.1v0$metab %in% output_sig_24$metab,]$log2FC
output_sig_24[output_sig_24$metab %in% df_BIT_24_10v0$metab,]$"24_10-0_log2FC" = df_BIT_24_10v0[df_BIT_24_10v0$metab %in% output_sig_24$metab,]$log2FC
output_sig_24[output_sig_24$metab %in% df_BIT_24_10v0.1$metab,]$"24_10-0.1_log2FC" = df_BIT_24_10v0.1[df_BIT_24_10v0.1$metab %in% output_sig_24$metab,]$log2FC

output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_24v0$metab,]$"untreated_24-0_log2FC" = df_BIT_untreated_24v0[df_BIT_untreated_24v0$metab %in% output_sig_untreated$metab,]$log2FC
output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_5v0$metab,]$"untreated_5-0_log2FC" = df_BIT_untreated_5v0[df_BIT_untreated_5v0$metab %in% output_sig_untreated$metab,]$log2FC
output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_5v24$metab,]$"untreated_5-24_log2FC" = df_BIT_untreated_5v24[df_BIT_untreated_5v24$metab %in% output_sig_untreated$metab,]$log2FC
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_24v0$metab,]$"persister_24-0_log2FC" = df_BIT_persister_24v0[df_BIT_persister_24v0$metab %in% output_sig_persister$metab,]$log2FC
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_5v0$metab,]$"persister_5-0_log2FC" = df_BIT_persister_5v0[df_BIT_persister_5v0$metab %in% output_sig_persister$metab,]$log2FC
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_5v24$metab,]$"persister_5-24_log2FC" = df_BIT_persister_5v24[df_BIT_persister_5v24$metab %in% output_sig_persister$metab,]$log2FC
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_24v0$metab,]$"dead_24-0_log2FC" = df_BIT_dead_24v0[df_BIT_dead_24v0$metab %in% output_sig_dead$metab,]$log2FC
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_5v0$metab,]$"dead_5-0_log2FC" = df_BIT_dead_5v0[df_BIT_dead_5v0$metab %in% output_sig_dead$metab,]$log2FC
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_5v24$metab,]$"dead_5-24_log2FC" = df_BIT_dead_5v24[df_BIT_dead_5v24$metab %in% output_sig_dead$metab,]$log2FC

output_sig_0[output_sig_0$metab %in% df_BIT_0_0.1v0$metab,]$"0_0.1-0_p" = df_BIT_0_0.1v0[df_BIT_0_0.1v0$metab %in% output_sig_0$metab,]$p
output_sig_0[output_sig_0$metab %in% df_BIT_0_10v0$metab,]$"0_10-0_p" = df_BIT_0_10v0[df_BIT_0_10v0$metab %in% output_sig_0$metab,]$p
output_sig_0[output_sig_0$metab %in% df_BIT_0_10v0.1$metab,]$"0_10-0.1_p" = df_BIT_0_10v0.1[df_BIT_0_10v0.1$metab %in% output_sig_0$metab,]$p
output_sig_5[output_sig_5$metab %in% df_BIT_5_0.1v0$metab,]$"5_0.1-0_p" = df_BIT_5_0.1v0[df_BIT_5_0.1v0$metab %in% output_sig_5$metab,]$p
output_sig_5[output_sig_5$metab %in% df_BIT_5_10v0$metab,]$"5_10-0_p" = df_BIT_5_10v0[df_BIT_5_10v0$metab %in% output_sig_5$metab,]$p
output_sig_5[output_sig_5$metab %in% df_BIT_5_10v0.1$metab,]$"5_10-0.1_p" = df_BIT_5_10v0.1[df_BIT_5_10v0.1$metab %in% output_sig_5$metab,]$p
output_sig_24[output_sig_24$metab %in% df_BIT_24_0.1v0$metab,]$"24_0.1-0_p" = df_BIT_24_0.1v0[df_BIT_24_0.1v0$metab %in% output_sig_24$metab,]$p
output_sig_24[output_sig_24$metab %in% df_BIT_24_10v0$metab,]$"24_10-0_p" = df_BIT_24_10v0[df_BIT_24_10v0$metab %in% output_sig_24$metab,]$p
output_sig_24[output_sig_24$metab %in% df_BIT_24_10v0.1$metab,]$"24_10-0.1_p" = df_BIT_24_10v0.1[df_BIT_24_10v0.1$metab %in% output_sig_24$metab,]$p

output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_24v0$metab,]$"untreated_24-0_p" = df_BIT_untreated_24v0[df_BIT_untreated_24v0$metab %in% output_sig_untreated$metab,]$p
output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_5v0$metab,]$"untreated_5-0_p" = df_BIT_untreated_5v0[df_BIT_untreated_5v0$metab %in% output_sig_untreated$metab,]$p
output_sig_untreated[output_sig_untreated$metab %in% df_BIT_untreated_5v24$metab,]$"untreated_5-24_p" = df_BIT_untreated_5v24[df_BIT_untreated_5v24$metab %in% output_sig_untreated$metab,]$p
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_24v0$metab,]$"persister_24-0_p" = df_BIT_persister_24v0[df_BIT_persister_24v0$metab %in% output_sig_persister$metab,]$p
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_5v0$metab,]$"persister_5-0_p" = df_BIT_persister_5v0[df_BIT_persister_5v0$metab %in% output_sig_persister$metab,]$p
output_sig_persister[output_sig_persister$metab %in% df_BIT_persister_5v24$metab,]$"persister_5-24_p" = df_BIT_persister_5v24[df_BIT_persister_5v24$metab %in% output_sig_persister$metab,]$p
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_24v0$metab,]$"dead_24-0_p" = df_BIT_dead_24v0[df_BIT_dead_24v0$metab %in% output_sig_dead$metab,]$p
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_5v0$metab,]$"dead_5-0_p" = df_BIT_dead_5v0[df_BIT_dead_5v0$metab %in% output_sig_dead$metab,]$p
output_sig_dead[output_sig_dead$metab %in% df_BIT_dead_5v24$metab,]$"dead_5-24_p" = df_BIT_dead_5v24[df_BIT_dead_5v24$metab %in% output_sig_dead$metab,]$p

write.csv(output_sig_0, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_welchANOVA_games_output_sig_0.csv'))
write.csv(output_sig_5, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_welchANOVA_games_output_sig_5.csv'))
write.csv(output_sig_24, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_welchANOVA_games_output_sig_24.csv'))
write.csv(output_sig_untreated, file.path(dir, 'analysis' ,'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_output_sig_untreated.csv'))
write.csv(output_sig_persister, file.path(dir, 'analysis' ,'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_output_sig_persister.csv'))
write.csv(output_sig_dead, file.path(dir, 'analysis', 'METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_output_sig_dead.csv'))

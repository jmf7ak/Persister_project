#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### set-up working environment #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

# set directory
setwd('..')
dir <- getwd()

# load the datasets
df_BIT_uniStats_0 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_0.csv"))
df_BIT_uniStats_5 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_5.csv"))
df_BIT_uniStats_24 <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_24.csv"))
df_BIT_uniStats_untreated <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_untreated.csv"))
df_BIT_uniStats_persister <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_persister.csv"))
df_BIT_uniStats_dead <- read.csv(file.path(dir, 'analysis', "METABOLON_DATA_origscaled_R_NoSampleNorm_greenhouseANOVA_games_df_BIT_uniStats_dead.csv"))

df_BIT_uniStats_0 <- select(df_BIT_uniStats_0, -X)
df_BIT_uniStats_5 <- select(df_BIT_uniStats_5, -X)
df_BIT_uniStats_24 <- select(df_BIT_uniStats_24, -X)
df_BIT_uniStats_untreated <- select(df_BIT_uniStats_untreated, -X)
df_BIT_uniStats_persister <- select(df_BIT_uniStats_persister, -X)
df_BIT_uniStats_dead <- select(df_BIT_uniStats_dead, -X)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining metabolites for 0-hour models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# there are not any significant differences between the untreated, persister and dead conditions at time 0
# these models should all have the same access/production of the same metabolites
# should be able to determine what these metabolites might be based on the individual condition specific comparisons evaluating differences between time

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining condition-specific metabolites for 5-hour models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## The purpose of this section is to determine which metabolites should be considered for integration

## metabs unique to untreated
untreated_unique_5_toP <- filter(df_BIT_uniStats_5, condition == '5_0.1-0' & p < 0.05 & log2FC < 0)
untreated_unique_5_toD <- filter(df_BIT_uniStats_5, condition == '5_10-0' & p < 0.05 & log2FC < 0)

length(intersect(untreated_unique_5_toP$metab, untreated_unique_5_toD$metab)) # 39 metabolites unique to untreated
untreated_unique_5 <- df_BIT_uniStats_5[df_BIT_uniStats_5$metab %in% intersect(untreated_unique_5_toP$metab, untreated_unique_5_toD$metab),]

## metabs unique to persister
persister_unique_5_toU <- filter(df_BIT_uniStats_5, condition == '5_0.1-0' & p < 0.05 & log2FC > 0)
persister_unique_5_toD <- filter(df_BIT_uniStats_5, condition == '5_10-0.1' & p < 0.05 & log2FC < 0)

length(intersect(persister_unique_5_toU$metab, persister_unique_5_toD$metab)) # 46 metabolites unique to persister
persister_unique_5 <- df_BIT_uniStats_5[df_BIT_uniStats_5$metab %in% intersect(persister_unique_5_toU$metab, persister_unique_5_toD$metab),]

## metabs unique to dead
dead_unique_5_toU <- filter(df_BIT_uniStats_5, condition == '5_10-0' & p < 0.05 & log2FC > 0)
dead_unique_5_toP <- filter(df_BIT_uniStats_5, condition == '5_10-0.1' & p < 0.05 & log2FC > 0)

length(intersect(dead_unique_5_toU$metab, dead_unique_5_toP$metab)) # 79 metabolites unique to dead
dead_unique_5 <- df_BIT_uniStats_5[df_BIT_uniStats_5$metab %in% intersect(dead_unique_5_toU$metab, dead_unique_5_toP$metab),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining condition-specific metabolites for 24-hour models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## The purpose of this section is to determine which metabolites should be considered for integration

## metabs unique to untreated
untreated_unique_24_toP <- filter(df_BIT_uniStats_24, condition == '24_0.1-0' & p < 0.05 & log2FC < 0)
untreated_unique_24_toD <- filter(df_BIT_uniStats_24, condition == '24_10-0' & p < 0.05 & log2FC < 0)

length(intersect(untreated_unique_24_toP$metab, untreated_unique_24_toD$metab)) # 113 metabolites unique to untreated
untreated_unique_24 <- df_BIT_uniStats_24[df_BIT_uniStats_24$metab %in% intersect(untreated_unique_24_toP$metab, untreated_unique_24_toD$metab),]

## metabs unique to persister
persister_unique_24_toU <- filter(df_BIT_uniStats_24, condition == '24_0.1-0' & p < 0.05 & log2FC > 0)
persister_unique_24_toD <- filter(df_BIT_uniStats_24, condition == '24_10-0.1' & p < 0.05 & log2FC < 0)

length(intersect(persister_unique_24_toU$metab, persister_unique_24_toD$metab)) # 40 metabolites unique to persister
persister_unique_24 <- df_BIT_uniStats_24[df_BIT_uniStats_24$metab %in% intersect(persister_unique_24_toU$metab, persister_unique_24_toD$metab),]

## metabs unique to dead
dead_unique_24_toU <- filter(df_BIT_uniStats_24, condition == '24_10-0' & p < 0.05 & log2FC > 0)
dead_unique_24_toP <- filter(df_BIT_uniStats_24, condition == '24_10-0.1' & p < 0.05 & log2FC > 0)

length(intersect(dead_unique_24_toU$metab, dead_unique_24_toP$metab)) # 49 metabolites unique to dead
dead_unique_24 <- df_BIT_uniStats_24[df_BIT_uniStats_24$metab %in% intersect(dead_unique_24_toU$metab, dead_unique_24_toP$metab),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### universal condition-specific metabolites #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untreated_unique_universal <- intersect(untreated_unique_5$metab, untreated_unique_24$metab) # 18
persister_unique_universal <- intersect(persister_unique_5$metab, persister_unique_24$metab) # 20
dead_unique_universal <- intersect(dead_unique_5$metab, dead_unique_24$metab) # 39

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining metabolites for untreated models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

untreated_unique_5_produced <- df_BIT_uniStats_untreated[df_BIT_uniStats_untreated$metab %in% untreated_unique_5$metab,] 
length(unique(untreated_unique_5_produced$metab))
untreated_unique_5_produced <- filter(untreated_unique_5_produced, condition == "untreated_5-0" & p < 0.05 & log2FC > 0)

untreated_unique_24_produced <- df_BIT_uniStats_untreated[df_BIT_uniStats_untreated$metab %in% untreated_unique_24$metab,] 
length(unique(untreated_unique_24_produced$metab))
untreated_unique_24_produced <- filter(untreated_unique_24_produced, condition == "untreated_5-24" & p < 0.05 & log2FC < 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining metabolites for persister models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

persister_unique_5_produced <- df_BIT_uniStats_persister[df_BIT_uniStats_persister$metab %in% persister_unique_5$metab,] 
length(unique(persister_unique_5_produced$metab))
persister_unique_5_produced <- filter(persister_unique_5_produced, condition == "persister_5-0" & p < 0.05 & log2FC > 0)

persister_unique_24_produced <- df_BIT_uniStats_persister[df_BIT_uniStats_persister$metab %in% persister_unique_24$metab,] 
length(unique(persister_unique_24_produced$metab))
persister_unique_24_produced <- filter(persister_unique_24_produced, condition == "persister_5-24" & p < 0.05 & log2FC < 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### determining metabolites for dead models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dead_unique_5_produced <- df_BIT_uniStats_dead[df_BIT_uniStats_dead$metab %in% dead_unique_5$metab,] 
length(unique(dead_unique_5_produced$metab))
dead_unique_5_produced <- filter(dead_unique_5_produced, condition == "dead_5-0" & p < 0.05 & log2FC > 0)

dead_unique_24_produced <- df_BIT_uniStats_dead[df_BIT_uniStats_dead$metab %in% dead_unique_24$metab,] 
length(unique(dead_unique_24_produced$metab))
dead_unique_24_produced <- filter(dead_unique_24_produced, condition == "dead_5-24" & p < 0.05 & log2FC < 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### output MetabsForIntegration #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(untreated_unique_5_produced, file.path(dir, 'analysis', 'MetabsForIntegration_untreated-unique-5-produced.csv'))
write.csv(untreated_unique_24_produced, file.path(dir, 'analysis', 'MetabsForIntegration_untreated-unique-24-produced.csv'))
write.csv(persister_unique_5_produced, file.path(dir, 'analysis', 'MetabsForIntegration_persister-unique-5-produced.csv'))
write.csv(persister_unique_24_produced, file.path(dir, 'analysis', 'MetabsForIntegration_persister-unique-24-produced.csv'))
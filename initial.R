# This will be my beginning of addressing this fun RNA-Seq challenge
# I have like a 20 step plan to go through, so here we go

# Call every library I can think of potentially using:


library(dplyr) 
library(data.table)
library(tidyr)
library(ggplot2) 
library(stringr)
library(knitr)
library(clipr) 

setwd("~/rnaseq")
rna_data <- read.csv("tcga.brca.rsem.csv",
                           stringsAsFactors = FALSE) #leave this is as false to start

#head(rna_data)

# 0. Check for NA -------

#any(is.na(rna_data))

# Returns FALSE, all fields have value present - low counts will be removed later

# 1. Filter for only primary tumor data source -------
#    can come back and implement design matrix
#    for other tumors in same patient later, but for now just want one thing per patient

#combined_results_df %>%
#  filter(grepl('surface area = ', V1)) -> combined_results_df

rna_data %>%
  filter(grepl('Primary Tumor', sample_type)) -> rna_data_primary_tumor_only
# 1212 rows to 1093


# 2. Drop barcode and sample type columns now that we've filtered the latter
#    Leaving patient_id unaltered; when we transpose and normalize we may then
#    use this in a design matrix or pull for modeling, scoop up other data from clinic

rna_data_genes_only <- rna_data_primary_tumor_only %>% 
  select(-c("bcr_patient_barcode","sample_type"))

# 3. Transpose this dataframe so we can begin our work! Patient_id/libraries as columns

#testrna <- as.data.frame(rna_data_genes_only[,1:9])

# from: https://stackoverflow.com/questions/7970179/transposing-a-dataframe-maintaining-the-first-column-as-heading
#tmydf = setNames(data.frame(t(testrna[,-1])), testrna[,1])

# 3 last lines DO work, so applying the same process

rna_data_t = setNames(data.frame(t(rna_data_genes_only[,-1])), rna_data_genes_only[,1])

rna_data_t_rd <- round(rna_data_t, digits = 0)

#tmydf_rounded <- round(tmydf, digits = 0)

#x <- as.data.frame(t(testrna))
#y <- as.data.frame()
#rna_data_t <- t(rna_data_genes_only)


#rna_data_only_genes <- subset(rna_data, select = -c("patient_id","bcr_patient_code","sample_type"))

# do we transpose yet? transpose is necessary for PCA
# No. Right now we need to get this normalized, so summing everything but the first 3 columns

# ok let's try dropping then bring back those columns

#patient_tumor_columns <- rna_data %>%
#  select(c("patient_id","bcr_patient_barcode","sample_type"))

#rna_data_only_genes <- rna_data %>% 
#  select(-c("patient_id","bcr_patient_barcode","sample_type"))

# Normalize ---------
# I'll adjust for read depth (total row counts), but regular normalization doesn't seem possible;
# I guess gene lengths could be imported to account for that, but this may have already been done during
# generation of the data. Find out how this was generated!!

#d <- t(apply(rna_data_only_genes, 1, function(x) x/sum(x)*1e6)) # this seems to work to divide by the sum of ea. row

# however, is this valuable??...
#


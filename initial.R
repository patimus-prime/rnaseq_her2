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
library(edgeR)  # NOTE: Must capitalize R lolololol


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

rna_data_primary_tumor_only <- rna_data %>%
  filter(grepl('Primary Tumor', sample_type)) 
# 1212 rows to 1093


# 2. Drop barcode and sample type columns now that we've filtered the latter ------------
#    Leaving patient_id unaltered; when we transpose and normalize we may then
#    use this in a design matrix or pull for modeling, scoop up other data from clinic

rna_data_genes_only <- rna_data_primary_tumor_only %>% 
  select(-c("bcr_patient_barcode","sample_type"))

# 3. Transpose this dataframe so we can begin our work! Patient_id/libraries as columns ----------------

#testrna <- as.data.frame(rna_data_genes_only[,1:9])

# from: https://stackoverflow.com/questions/7970179/transposing-a-dataframe-maintaining-the-first-column-as-heading
#tmydf = setNames(data.frame(t(testrna[,-1])), testrna[,1])

# 3 last lines DO work, so applying the same process

rna_data_t = setNames(data.frame(t(rna_data_genes_only[,-1])), rna_data_genes_only[,1])


# 4. Round up all data ----------------

rna_data_t_rd <- round(rna_data_t, digits = 0)
# hold for data from other files
# this will be principally used for most stuff in edgeR, but
# we'll want to pivot back AFTER normalizing and use HER2 column and others for modeling


# 5. We must get data from other sources --------------------
#    so we can actually verify that the 
#    way patients' expression differs is actually due to +/- 

#    A. Get patient ID and IHC, remove NA; we'll use this to define groups
#       We'll use copy number later when modeling, more useful for verifying

      clinical_data <- read.csv("brca_tcga_clinical_data.csv",
                                stringsAsFactors = FALSE)
      
      patient_IHC <- clinical_data %>%
        select(Patient.ID, IHC.HER2)
      # WARNING: contains NA, Intermediate, Equivocal - perhaps useful in a design matrix later.
      # maybe we can also come back and grab copy number to see how that clusters!
      # will require definition of design matrix

#     B. Grab the positive and negative groups.
      
      positive_patients <- patient_IHC %>%
        filter(grepl('Positive', IHC.HER2)) 
      #164 patients
      negative_patients <- patient_IHC %>%
        filter(grepl('Negative', IHC.HER2))
      #567 patients

# now we have our groups defined, can specify for edgeR below

# 6. Into edgeR for normalization
      
      
      y <- drop_na(patient_IHC)
      
      matching_ids <- intersect(y$Patient.ID,colnames(rna_data_t_rd,))
      # therefore mathcing records contains all paitent ids with both IHC results and RNASeqs of primary tumors
      
      #matching_ids <- as.data.frame(matching_ids) 
      
       
      
      # so let's dump what's not in here for our below analysis
      
      xx <- select(rna_data_t_rd, matching_ids)
      
      # now a discrepancy between rna_data_t_rd and y, the no NA patient IHC data)
      # so like... I guess 
      
      #u <- as.data.frame(matching_ids)
    
    Patient.ID <- intersect(y$Patient.ID, matching_ids) # NEW MATCHING IDS
    newMATCHES <- as.data.frame(Patient.ID)
    #final <- as.data.frame(uu)
    
   # only keeps the values in newMATCHES, left df, if extras detected
    EXALTATION <- y %>%
      semi_join(newMATCHES, by = "Patient.ID")
    
    weird_rows <- y %>%
      anti_join(newMATCHES, by = "Patient.ID") # FIXME: need to investigate this later
    
    EXCELLENCE <- EXALTATION %>%
      anti_join(weird_rows, by = "Patient.ID")
    
    
    # YES!!!!!!!
    eee <- distinct(EXALTATION) #weird lol
    
    # ok so we now have all columns that match, after much derision and confusion
      
      
      # now let's get the ihc that matches up
      d <- DGEList(counts = xx, group = factor(eee$IHC.HER2))
      d
      dim(d) #cool, everything is working
      
      
      d <- calcNormFactors(d)
      
      d$samples
      
      ################
      # BEHOLD! THE REMAINDER IS SIMPLE NOW THIS OBJECT IS DEFINED. WE MAY NOW GET CLUSTERS ETC.
      # THEN WE CAN TAKE THE OUTPUT NORMALIZED FACTORS, THE MOST CRITICAL COMPONENTS AS DETERMINED 
      # BY PCA, AND THEN FEED THAT INTO THE LOGISTIC REGRESSION MODEL ALONG WITH THE ALLELE/COUNT THING
      #############
      
      
      # cool
      # BOOM THAT'S FINALLY DEFINED WHAT WHAT!!!!!
      
      #dge <- calcNormFactors(dge, method = "TMM")
      #tmm <- cpm(dge)
      
      norm_counts <- cpm(d) # interesting how much stuff changed lol, but those were raw counts and these are cpm from TMM
      
      
      
      # INCREDIBLE -----------------------------
      
      # these are the rnaseq results that match up with the shit
      
      # filter the rows of genes that have very low expression 
      # use cutoff of 1 counts per million; i'll try a few things but 
      # a quick google wasn't enough to find a good rule of thumb, so going off edgeR guide
      cutoff = 1 # counts per million
      keptGenes <- rowSums(cpm(rna_data_t_rd) > cutoff) >= ncol(rna_data_t_rd) #LAST PART = # LIBRARIES, SO 
                                                                            # THE ENTIRE ROW HAS AT LEAST 1 READ/CPM PER LIBRARY
      
      # now i think I might be skipping some steps... let's eat first 
      # ok, let's get this stuff into a list and make sure we're starting from the beginning
      super <- rna_data_t_rd[keptGenes,]
      super <- calcNormFactors(super, method = "upperquartile")
      
     
# Normalize ---------
# I'll adjust for read depth (total row counts), but regular normalization doesn't seem possible;
# I guess gene lengths could be imported to account for that, but this may have already been done during
# generation of the data. Find out how this was generated!!

#d <- t(apply(rna_data_only_genes, 1, function(x) x/sum(x)*1e6)) # this seems to work to divide by the sum of ea. row

# however, is this valuable??...
#


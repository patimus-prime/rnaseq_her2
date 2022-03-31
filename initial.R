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
library(tibble)
library(magrittr)

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
      # Also this turns out to be useless, but maybe can use later
      
      positive_patients <- patient_IHC %>%
        filter(grepl('Positive', IHC.HER2)) 
      #164 patients
      negative_patients <- patient_IHC %>%
        filter(grepl('Negative', IHC.HER2))
      #567 patients

#     C. Not sure we got the groups correct, need to specify and 
      # match up everything to reference the way edgeR does it.
      # BEWARE: THIS IS CONFUSING BUT WORKS. LOOK INTO ALTERNATIVE LATER.
      # PROPOSAL: USE DISTINCT() UPSTREAM, THEN PERHAPS WE CAN AVOID THE BACK AND FORTH
      #           ... MERGE IS LIKELY NECESSARY THOUGH.
      
      # get only the patients that have IHC results
      y <- drop_na(patient_IHC)
      matching_ids <- intersect(y$Patient.ID,colnames(rna_data_t_rd,))
      # therefore mathcing records contains all paitent ids with both IHC results and RNASeqs of primary tumors
    
      # get only the rna data that is from patients with IHC results
      matched_rna_data <- select(rna_data_t_rd, matching_ids)
      
      # but awkwardly, we now have to use the one IDs that are shared
      # to go back to our IDs and IHC results to see what matches to RNA.
      # This seems excessively backwards and I'll take another look if time permits
    
      Patient.ID <- intersect(y$Patient.ID, matching_ids) # NEW MATCHING IDS
      newMATCHES <- as.data.frame(Patient.ID)
    
   # only keeps the values in newMATCHES, left df, if extras detected
    matched_patient_IHC <- y %>%
      semi_join(newMATCHES, by = "Patient.ID")
    
    weird_rows <- y %>%
      anti_join(newMATCHES, by = "Patient.ID") # FIXME: need to investigate this later
    
    # Remove duplicates -- look into the weird_rows above if time permits
    matched_patient_IHC <- distinct(matched_patient_IHC) #weird lol
    
    # EXALT!! so we now have all columns that match, after much derision and confusion

# 6. Into edgeR for normalization
    
    dim(matched_rna_data) # 20502 genes
    
    # A. First filter out genes expressed very little
    
    cutoff = 1 # counts per million; other cutoffs are used but this is what we'll start with
    
    keptGenes <- rowSums(cpm(matched_rna_data) > cutoff) >= ncol(matched_rna_data) #LAST PART == # LIBRARIES, SO 
    #THE ENTIRE ROW HAS AT LEAST 1 READ/CPM PER LIBRARY
    
    matched_rna_data <- matched_rna_data[keptGenes,]
    
    dim(matched_rna_data) # 9105 genes with not very little expression after filtering
                          # in case of confusion: we are removing columns (genes, counts), not patients
    
    # B. Create DGE List used in all of edgeR functions 
    
      # define the edgeR list object that matches everything up and makes life easy
      d <- DGEList(counts = matched_rna_data, group = factor(matched_patient_IHC$IHC.HER2))
      d
      dim(d) #cool, everything is working
      
    # C. Normalization!!!!! CRITICAL STEP! ###########
      # perform normalization, get factors based on edgeR's normalization algorithm (excluding VERY lowly expressed genes)
      # btw edgeR uses TMM and a bunch of funky fresh math to do this
      d <- calcNormFactors(d)
      
      d$samples # more checking everything's A-OK
      
      ################
      # BEHOLD! THE REMAINDER IS SIMPLE NOW THIS OBJECT IS DEFINED, NORMALIZED. WE MAY NOW GET CLUSTERS ETC.
      # THEN WE CAN TAKE THE OUTPUT NORMALIZED FACTORS, THE MOST CRITICAL COMPONENTS AS DETERMINED 
      # BY PCA, AND THEN FEED THAT INTO THE LOGISTIC REGRESSION MODEL ALONG WITH THE ALLELE/COUNT THING
      #############
      
   # D. Get normalized counts, can filter for most critical genes later and input to logistic regression model
      norm_counts <- cpm(d) # Data look different,but those were raw counts and these are cpm from TMM
      
      # from: https://www.biostars.org/p/317701/
      #dge <- calcNormFactors(dge, method = "TMM")
      #tmm <- cpm(dge)
      
      
# 7. Further processes in edgeR ----------------
      # without logs, very poor graph
      #plotMDS(matched_rna_data)
      
      # with logs, use only 500 top genes, use group labels instead of patient IDS
      #plotMDS.DGEList(matched_rna_data, top = 500, labels = d$samples$group, gene.selection = "common", pch = ".") #,  pch = points(d$samples$group))
      
      mds <- plotMDS.DGEList(matched_rna_data, top = 500, 
                           method = "logFC"  , labels = d$samples$group, 
                           gene.selection = "common", var.explained = TRUE)
    
      
      uberplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(d$samples$group))
      ggplot(uberplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      
     #  mds <- plotMDS(x)
     #  toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(paste0("Grp", rep(1:2, each = 3))))
     #  library(ggplot2)
     #  ggplot(toplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
     #  
     #  toplot <- as.data.frame(mds@.Data[[3]]) %>%
     #    left_join(d$samples$group)
     # #mds@.Data[[3]]
     #  mds@.Data[[3]] %>%
     #    as.data.frame() %>%
     #    #set_colnames(c("Dim1", "Dim2")) %>%
     #    #rownames_to_column("SampleID") %>%
     #    left_join(d$samples$group)
     #  
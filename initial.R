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
      
      
# 7. Attempt clustering in edgeR ----------------
      # without logs, very poor graph
      #plotMDS(matched_rna_data)
      
      # with logs, use only 500 top genes, use group labels instead of patient IDS
      #plotMDS.DGEList(matched_rna_data, top = 500, labels = d$samples$group, gene.selection = "common", pch = ".") #,  pch = points(d$samples$group))
      
      mds <- plotMDS.DGEList(d, top = 500, 
                           method = "logFC"  , labels = d$samples$group, 
                           gene.selection = "common", var.explained = TRUE)
    
      # comment these lines and run one at a time or R crashes idk man
      #uberplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(d$samples$group))
      #ggplot(uberplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      
      
# 8. Try other stuff -------
    # let's try differntially expressed genes - maybe there's just a bunch that confuddle things
      
      # d <- estimateCommonDisp(d, vebose = TRUE)
      # d <- estimateTagwiseDisp(d, trend = "none")
      # plotBCV(d, cex = 0.4)
      # et <- exactTest(d)
      # topTags(et, n = 20)
      # 
      # detags <- rownames(topTags(et, n=20))
      # df_detags <- cpm(d)[detags,]
      # df_detags <- as.data.frame(df_detags)
      
      # STILL NO PICKLES!!!!!!!!!!!!!!!
      
# 9. Ok, let's try specifying a comparison against positive and negative testers
      # et2 <- exactTest(d, pair = c("Positive","Negative"))
      # topTags(et2, n = 20)
      # # AHHHHHHH NOW WE FIND ERBB2 WHAT WHAT!!@!!!!!! AHAHAHAAHAHAHAHHAHA FUCK YES
      # 
      # # well, we can do a quick smear plot, but more interesting will be to go back
      # # and redo PCA, only going for the positive and negative
      # # this does suggest that this whole method can only be supplementary to ICH etc.
      # summary(de <- decideTestsDGE(et2, p=0.05))
      # detags <- rownames(d)[as.logical(de)]
      # plotSmear(et2, de.tags=detags)
      # abline(h = c(-2, 2), col = "blue")
      
      #  cool. but better if we could mark indiv. genes
      # that'll be a tomorrow problem to put in tags for ERBB2 etc.
      # for now let's go back do clustering and only do positive/negative
    
      
      # now THIS is really interesting. We don't even see HER2 on here. This suggests
      # maybe we have to go back and ONLY look at positive/negative relations.
      # equivocal etc. may affect the calculation of variance etc.
      
      
# 10. Only Positive/negative results into PCA. May do design matrix later
      # but for now this is the most interesting to go after. Can demonstrate
      # the model first before doing design matrix stuff - this would simply
      # be enumerating the equivocal etc. and therefore not be of additional use
      
      ####################3
      # A. Prepare the data for going into object
      ######################
      
      # this is going to require constructing a new dge object ... hopefully more straightforward this time
      # let's try appending the positive/negative dataframes first
      io_patients <- rbind(positive_patients, negative_patients)
      io_patients <- distinct(io_patients)
  
      io_rna_match <- intersect(io_patients$Patient.ID, colnames(rna_data_t_rd))
      # mismatch again... just use same method to eliminate the weird rows
      
      io_matched_rna <- select(rna_data_t_rd, io_rna_match)
    
      Patient.ID <- intersect(io_patients$Patient.ID, io_rna_match) # NEW MATCHING IDS
      
      ioMATCHES <- as.data.frame(Patient.ID)
      
      io_group <- io_patients %>%
        semi_join(ioMATCHES, by = "Patient.ID")
      # cool, easy this time.
      
      # so we have: 
      # io_group
      # and
      # io_matched_rna
      
      # now repeat the process same as before.
      
      ######################3
      # B. Create the object
      ###################3#####
      
      cutoff = 1 # counts per million; other cutoffs are used but this is what we'll start with
      
      keep <- rowSums(cpm(io_matched_rna) > cutoff) >= ncol(io_matched_rna) #LAST PART == # LIBRARIES, SO 
      #THE ENTIRE ROW HAS AT LEAST 1 READ/CPM PER LIBRARY
      
      io_matched_rna <- io_matched_rna[keep,]
      # down to 9289 genes
    
      io_dge <- DGEList(counts = io_matched_rna, group = factor(io_group$IHC.HER2))
      
      io_dge <- calcNormFactors(io_dge)
      
      io_dge$samples # more checking everything's A-OK
      
      io_norm_counts <- cpm(io_dge)
      
      io_mds <- plotMDS.DGEList(io_dge, top = 500, 
                                method = "logFC"  , labels = io_dge$samples$group, 
                                gene.selection = "common", var.explained = TRUE)
      
      
      # comment these lines and run one at a time or R crashes idk man
      ioplot <- data.frame(Dim1 = io_mds$x, Dim2 = io_mds$y, Group = factor(io_dge$samples$group))
      
      #ggplot(ioplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      
      # mds_bcv <- plotMDS.DGEList(d, top = 500, 
      #                        method = "bcv"  , labels = d$samples$group, 
      #                        gene.selection = "common", var.explained = TRUE)
      # bcvplot <- data.frame(Dim1 = mds_bcv$x, Dim2 = mds_bcv$y, Group = factor(d$samples$group))
      # ggplot(bcvplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      # 
      
      
      # 
      # 
      # mds_bcv <- plotMDS.DGEList(matched_rna_data, top = 500, 
      #                        method = "logFC"  , labels = d$samples$group, 
      #                        gene.selection = "common", var.explained = TRUE)
      # 
      # uberplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(d$samples$group))
      # ggplot(uberplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      # 
      # so basically we do get the plot we want, but... it is not very informative!!
      # RIP!
      
      # so now we move on to other analyses
      
# 8. WTF am i doing
      
      
      # mds$var.explained
      # mdsDf <- as.data.frame(mds@.Data[[3]])
      # left_join(mdsDf,as.data.frame(d$samples$group), by= character())
      # d$samples$group
      # mds
      # mdsDf$new <- d$samples$group
      # 
      # mds$eigen.vectors
      # mds$gene.selection  
      # mds$dim.plot
      # mds$
      
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
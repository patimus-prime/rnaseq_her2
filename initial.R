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
library(tibble)

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
    
      uberplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(d$samples$group))
      ggplot(uberplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      
      
# 8. Try other stuff -------
    # let's try differntially expressed genes - maybe there's just a bunch that confuddle things
      
      d <- estimateCommonDisp(d, vebose = TRUE)
      d <- estimateTagwiseDisp(d, trend = "none")
      # plotBCV(d, cex = 0.4)
      # et <- exactTest(d)
      # topTags(et, n = 20)
      # 
      # detags <- rownames(topTags(et, n=20))
      # df_detags <- cpm(d)[detags,]
      # df_detags <- as.data.frame(df_detags)
      
      # STILL NO PICKLES!!!!!!!!!!!!!!!
      
# 9. Ok, let's try specifying a comparison against positive and negative testers
      et2 <- exactTest(d, pair = c("Positive","Negative"))
      io_DEGdf <- topTags(et2, n = 10)
      # THERE'S ERBB2 LET'S GET A LIST OF THOSE GENES HOMIE!
      
      io_DEGdf <- as.data.frame(io_DEGdf)
      io_DEG <- tibble::rownames_to_column(io_DEGdf, "Genes")
      io_ls_DEG <- io_DEG$Genes
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
      
      ggplot(ioplot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      
      # this is more informative than the last - although not too much
      # there's a cluster to the right, but it appears to be a small proportion
      # of the overall data...
      
      # maybe if we only use ERRB2 data to separate things?
      # we can just plot based on ERRB2 expression. 
      # i'll do that, but otherwise move onto modeling
      
# 11. Let's just get the normalized expression data of the different thresholds
      # and see how it may differ. Simple barplot. Nothing fancy.
      
      # So...
      io_t <- as.data.frame(t(as.data.frame(io_norm_counts)))
      # this one here, is what you'd wanna filter down for only the ERBB2 gene
      # and cross-reference with the groups. dunno if possible outside the object
      # but anyway i'm fuckin done lol for tonight
      
      # from: https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
      # df <- tibble::rownames_to_column(df, "VALUE")
      
      #cool easy
      io_t_norm <- tibble::rownames_to_column(io_t, "Patient.ID")
      
      # now drop everything except  ERBB2 and Patient.ID
      # FIXME: UNNCESSARY PROCESS
      her2df <- io_t_norm[c("Patient.ID", "ERBB2")]
      
      # merge in the patient ID, not a DGE so require grouping info
      her2df <- her2df %>%
        left_join(io_group, by = "Patient.ID")
      
      #her2_positive <- her2df[her2df$IHC.HER2 == 'Positive',] # beware the comma!
      
      #her2_negative <- her2df[her2df$IHC.HER2 == 'Negative',]
      
      #ggplot(her2df, aes(x=as.factor(IHC.HER2), y = ERBB2)) +
      #         geom_bar(stat = "summary", position = "dodge") # stat = identity?
      
      # VERY INTERESTING!!!
      ggplot(her2df, aes(x = ERBB2, colour = IHC.HER2)) +
        geom_density()
      
      # this is so interesting in fact, we'll go back and do another plot like this
      # for all the IHC levels, and see how that does. 
      # then if we still see it, we can see about also plotting other DE genes for either df
      
      d_t <- as.data.frame(t(as.data.frame(norm_counts)))
      d_t <- tibble::rownames_to_column(d_t, "Patient.ID")
      # actually let's merge beforehand
      d_t <- d_t %>%
        left_join(matched_patient_IHC, by = "Patient.ID")
      
      ggplot(d_t, aes(x = ERBB2, colour = IHC.HER2)) +
        geom_density()
      
      
      io_ls_DEG
      
      ggplot(d_t, aes(x = STARD3, colour = IHC.HER2)) +
          geom_density()
      ggplot(d_t, aes(x = PGAP3, colour = IHC.HER2)) +
        geom_density()
      ggplot(d_t, aes(x = C17orf37, colour = IHC.HER2)) +
        geom_density()
      ggplot(d_t, aes(x = ORMDL3, colour = IHC.HER2)) +
        geom_density()
      
      
      # ok, let's try clustering with just the top 5, it doesn't look promising , but maybe
      # if not, we'll move onto the logit function and call it a day lol\
      
# 12. Reclustering with only the top 5 most DE genes and hoping for a miracle!
      
      # we'll do with +/- then all and see, 'cuz right now it's not looking too hot for using this as
      # the sole diagnostic tool
      #her2df <- io_t_norm[c("Patient.ID", "ERBB2")]
      top5df_io <- io_t_norm[c("Patient.ID","ERBB2","STARD3","PGAP3","C17orf37","ORMDL3")]
      top5df_io <- top5df_io %>%
        remove_rownames %>% column_to_rownames(var = "Patient.ID")
      top5df_io <- t(top5df_io)
      io2_dge <- DGEList(counts = top5df_io,group=factor(io_group$IHC.HER2))
      io2_dge <- calcNormFactors(io2_dge)
      
      io2_mds <- plotMDS.DGEList(io2_dge,
                                method = "logFC"  , labels = io2_dge$samples$group, 
                                gene.selection = "pairwise", var.explained = TRUE)
      
      # comment these lines and run one at a time or R crashes idk man
      io2plot <- data.frame(Dim1 = io2_mds$x, Dim2 = io2_mds$y, Group = factor(io2_dge$samples$group))
      
      ggplot(io2plot, aes(Dim1, Dim2, colour = Group)) + geom_point()
      # still not able to cluster these things apart, even with just using most
      # differentially expressed genes. rip.
      
# so let's try instead integrating as much as we can; let's put a csv together to feed into scikit learn

# 13. Outputting a CSV to use in scikit learn in python, just throwing as many features in there 
      # as I can think may help and which don't have NA
      # maybe use FISH here if it has enough rows
      # just gonna drop all the NA anyway
      io_csv_df <- as.data.frame(t(top5df_io))
      io_csv_df <- tibble::rownames_to_column(io_csv_df, "Patient.ID")
      
      io_csv_df <- io_csv_df %>%
        left_join(io_group, by = "Patient.ID")
      io_csv_df <- distinct(io_csv_df)
      
      # ok that's all gucci. get all data in and put io at end w dplyr
      
      # other data to grab:
      # Ethnicity
      # Fraction Genome Altered
      # Race Category
      # Sex
      # Age
      
      # method: get a df of these from clinical data
      # then merge via PatientID
      
      # a lot of NA in Race and Ethnicity, don't wanna bias so I'll drop those
      # and sex, which is majority women. So we're mostly using age, fraction genome 
      # and RNASeq data to predict IHC, not great, not terrible
      tocsv_clinical <- clinical_data[c("Patient.ID","Diagnosis.Age","Fraction.Genome.Altered")]
      
      all_to_csv <- io_csv_df %>%
        left_join(tocsv_clinical, by = "Patient.ID")
      
      all_to_csv <- distinct(all_to_csv)
      all_to_csv <- drop_na(all_to_csv)
      
      all_to_csv <- all_to_csv %>% 
        relocate(IHC.HER2, .after = last_col())
      # OK, now we can save to a csv, input to a python script, get the
      # model and somehow include in Rmd and write things up. WOO!
      write.csv(all_to_csv,"~/rnaseq/tologit.csv",row.names = FALSE) # leave pos/neg to be sure when building model
      rm(clinical_data)
      rm(temp)      
            
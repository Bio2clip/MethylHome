#!/usr/bin/env Rscript


# ----------------------------------------------------------------------
# Script : extract_other_qc_metrics.R
# Object  : Extract all QC
# Librairies : ewastools, dplyr 
# This script has been written using the script developped by Yvan Nicaise & Clementine Decamps, 2025-2026
# ----------------------------------------------------------------------

#### --- DOCUMENTATION --- ####
#' Compute Quality Control metrics
#' 
#' @Description 
#' 
#' Computes QC metrics as described in the 'BeadArray Controls Reporter Software Guide' from Illumina. 
#' Computes also log2Meth and log2 unmeth Intensities and Detection rate.
#' 
#' @Usage 
#' 
#' extract_qc_metrics(meth_rds, sample_id)
#' 
#' @Arguments 
#' 
#' meth_rds      Output of calling ewastools::read_idats. You can give raw or dye-biaised corrected data
#' sample_id     Name of the sample
#' 
#' @Value 
#' 
#' Return a dataframe containing all computed metrics per sample.
#' 
#' @Details
#' 
#' The metrics computed by the function are listed below : 
#' 
#'  - log2 Methylated intensity
#'  - log2 Unmethylated intensity
#'  - Detection Rate
#'  - Restoration
#'  - Staining Green
#'  - Staining Red
#'  - Extension Green
#'  - Extension Red
#'  - Hybridization High/Med
#'  - Hybridization Med/Low
#'  - Target Removal I
#'  - Target Removal II
#'  - Bisulfite Conversion I Green
#'  - Bisulfite Conversion I Green Bkg
#'  - Bisulfite Conversion I Red
#'  - Bisulfite Conversion I Red Bkg
#'  - Bisulfite Conversion II
#'  - Bisulfite Conversion II Bkg
#'  - Specificity I Green
#'  - Specificity I Red
#'  - Specificity II
#'  - Specificity II Bkg
#'  - Nonpolymorphic Green
#'  - Nonpolymorphic Red- Nonpolymorphic Red
#' 
#' @Note 
#' 
#' If we want to compute metrics as performed by 'BeadArray Controls Reporter Software' from Illumina
#' you have to give raw data before any normalization step to the function. 

#### --- CODE --- ####
library(stringr)
library(ewastools)
library(dplyr)
library(tidyr)

meth_QC <- readRDS("${meth_rds}")
sample_name <- "${sample_id}"

### --- QC metrics

# Throw an error if bad quality samples

QC <- control_metrics(meth_QC)
table_res <- c()
for(i in 1:length(QC)){
  loc_metrique = QC[[i]]
  table_res = rbind(table_res, c(names(QC)[i], attr(loc_metrique,"threshold"), loc_metrique))
    
}
colnames(table_res) <- c("Metrique", "Seuil", meth_QC[["meta"]][["sample_id"]])
 
table_res <- as.data.frame(table_res)
  
table_res_t <- t(table_res) # Transpose table to access more easily Metrics values
colnames(table_res_t) <- gsub(" ", "_", x=table_res[["Metrique"]]) # To prevent bad naming, e.g. using spaces
colnames(table_res_t) <- gsub('\\\\(Bkg)', "Bkg.", x = colnames(table_res_t)) # Better name for plotting (doesn't like "()")
table_res_t <- as.data.frame(table_res_t[-1,]) # get rid of unwanted values
table_res_t <- as.data.frame(table_res_t[-1,])
table_res_t[["Sample"]] <- rownames(table_res_t)
  
### --- Main QC Metrics 
  
# Log2 (Un)Meth Intensities
Log2MethIntensity_ewas   <- round(log2(matrixStats::colMedians(meth_QC[["M"]]  + 1, na.rm=TRUE)), 3)
Log2UnmethIntensity_ewas   <- round(log2(matrixStats::colMedians(meth_QC[["U"]]  + 1, na.rm=TRUE)), 3)
  
names(Log2MethIntensity_ewas) <- meth_QC[["meta"]][["sample_id"]]  # Name the intensities according to their sample of origin
names(Log2UnmethIntensity_ewas) <- meth_QC[["meta"]][["sample_id"]]

# Check how many probes are completely missing (M + U)
if (colSums(is.na(meth_QC[["M"]] + meth_QC[["U"]])) > 1000) {
  
  # DetectionP computation
  detP_ewas <- NA
  
  # Create summary dataframe
  main_qc_df <- data.frame(
    Sample_Name = sample_name,
    Sample = names(Log2MethIntensity_ewas),
    DetectionRate = NA,
    Log2MethIntensity = Log2MethIntensity_ewas,
    Log2UnmethIntensity = Log2UnmethIntensity_ewas,
    check.names = FALSE
  ) |>
    dplyr::left_join(table_res_t, by = "Sample")
  
} else{
  # DetectionP computation
  detP_ewas <- meth_QC %>% ewastools::detectionP()
  
  # Create summary dataframe
  main_qc_df <- data.frame(
    Sample_Name = sample_name,
    Sample = names(Log2MethIntensity_ewas),
    DetectionRate = round(colMeans(detP_ewas[["detP"]] < 0.05, na.rm=TRUE), 4),
    Log2MethIntensity = Log2MethIntensity_ewas,
    Log2UnmethIntensity = Log2UnmethIntensity_ewas,
    check.names = FALSE
  ) |>
    dplyr::left_join(table_res_t, by = "Sample")
  
}

### -- Export all QC metrics in a csv file
write.table(main_qc_df, paste0(sample_name, "_qc_metrics_output.tsv"), row.names = F, sep = "\t", quote = FALSE, dec = ".")

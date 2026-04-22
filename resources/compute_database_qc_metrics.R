#!/usr/bin/env Rscript


# ----------------------------------------------------------------------
# Script : compute_database_qc_metrics.R
# Object  : Extract and store all QC from a reference dataset
# Librairies : ewastools, stringr, dplyr, tidyr
# This script has been written using the script developped by Yvan Nicaise & Clementine Decamps, 2025-2026
# ----------------------------------------------------------------------


#### --- DOCUMENTATION --- ####
#' Compute and write Quality Control metrics for given IDAT files
#' 
#' @Description 
#' 
#' Computes QC metrics as described in the 'BeadArray Controls Reporter Software Guide' from Illumina. 
#' Computes also log2Meth and log2 unmeth Intensities and Detection rate.
#' This script can be used to create a reference file to launch the pipeline upon. 
#' 
#' @Usage 
#' 
#' Rscript compute_database_qc_metrics.R idat_qc_folder output_directory
#' 
#' @Arguments 
#' 
#' idat_dir      directory of the idat files from the reference database
#' outdir        directory on which to write tsv output file
#' 
#' @Value 
#' 
#' Write a tsv file containing all computed metrics per sample.
#' 
#' @Details
#' 
#' The metrics computed by the function are listed below : 
#' 
#'  - Detection Rate
#'  - log2 Methylated intensity
#'  - log2 Unmethylated intensity
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

### --- Read idat files 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("stop Usage : Rscript compute_database_qc_metrics.R <idat_directory> <output_directory>")
}

idat_dir <- normalizePath(args[1])
outdir <- args[2]


# Get rid of "_Grn.idat" 

files = list.files(idat_dir, full.names = T)
files = files[stringr::str_detect(files, "\\.idat")]
files = unique(stringr::str_remove_all(files, "_Red\\.idat|_Grn\\.idat"))

print(paste("Sample prefix: ", files))

meth_QC <- read_idats(files, quiet=FALSE) # `quiet=TRUE` suppresses the progress bar

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
colnames(table_res_t) <- gsub('\\(Bkg)', "Bkg.", x = colnames(table_res_t)) # Better name for plotting (doesn't like "()")
table_res_t <- as.data.frame(table_res_t[-1,]) # get rid of unwanted values
table_res_t <- as.data.frame(table_res_t[-1,])
table_res_t[["Sample"]] <- rownames(table_res_t)

### --- Main QC Metrics 

# Log2 (Un)Meth Intensities
Log2MethIntensity_ewas   <- round(log2(matrixStats::colMedians(meth_QC[["M"]]  + 1, na.rm=TRUE)), 3)
Log2UnmethIntensity_ewas   <- round(log2(matrixStats::colMedians(meth_QC[["U"]]  + 1, na.rm=TRUE)), 3)

names(Log2MethIntensity_ewas) <- meth_QC[["meta"]][["sample_id"]]  # Name the intensities according to their sample of origin
names(Log2UnmethIntensity_ewas) <- meth_QC[["meta"]][["sample_id"]]

# DetectionP computation
detP_ewas <- meth_QC %>% ewastools::detectionP.neg()

# Create summary dataframe
main_qc_df <- data.frame(
  Sample = names(Log2MethIntensity_ewas),
  DetectionRate = round(colMeans(detP_ewas[["detP"]] < log10(0.05), na.rm=TRUE), 4),
  Log2MethIntensity = Log2MethIntensity_ewas,
  Log2UnmethIntensity = Log2UnmethIntensity_ewas,
  check.names = FALSE
) |>
  dplyr::left_join(table_res_t, by = "Sample")

### -- Export all QC metrics in a csv file
write.table(main_qc_df, paste0(outdir, "/qc_metrics_output_db.tsv"), row.names = F, sep = "\t", quote = FALSE, dec = ".")

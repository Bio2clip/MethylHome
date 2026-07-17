#!/usr/bin/env Rscript


# Compute xy methylation intensities
# Libraries : ewastools, stringr, forcats, dplyr, tidyr

#### --- DOCUMENTATION --- ####
#' Compute X and Y intensities from methylation data
#' 
#' @Description 
#' 
#'Compute X and Y intensities from methylation data and store them in a tsv file.
#'The input file is based on BeadArray sample_sheet with columns Sample_IDAT and Gender
#' 
#' @Usage 
#' 
#' compute_xy_intensities.R(sample_info_df, meth_rds, sample_id)
#' 
#' @Arguments 
#' 
#' sex                     Gender of the studied sample
#' meth_QC                 R data object containing a meth object as provided by ewastools::read_idat function
#' sample_id               name of the given samples
#' 
#' 
#### --------------------------- ####


library(forcats)
library(dplyr)
library(tidyr)
library(stringr)
library(ewastools)


sex <- "${sex}"

meth_QC <- readRDS("${meth_rds}")
sample_name <- "${sample_id}"

# Retrieve threshold for bad quality samples
threshold <- ${quality_threshold}

### --- Put capital letter to lower letter
if (sex == "M" || sex == "Male") {
  sex = "m"
} else if (sex == "F" || sex == "Female") {
  sex = "f"
}

#### --- Compute X and Y intensities
if (colSums(is.na(meth_QC[["M"]] + meth_QC[["U"]])) <= threshold){
  
  # Apply dye-biaised correction
  meth_extra_QC_corrected <- meth_QC %>% correct_dye_bias 
  # Extract normalized average X and Y intensities
  sex_info_xy <- ewastools::check_sex(meth_extra_QC_corrected)
  
  # Write df 
  sex_info <- data.frame(Sample_Name = sample_name,
                         Sample_IDAT = meth_QC[["meta"]][["sample_id"]],
                         X = sex_info_xy[["X"]],
                         Y = sex_info_xy[["Y"]],
                         Gender = sex)
  
} else {
  sex_info <- data.frame(
    Sample_Name = sample_name,
    Sample_IDAT = meth_QC[["meta"]][["sample_id"]],
    X = NA,
    Y = NA, 
    Gender = "U"
  )
}

# Save computed X and Y intensities
write.table(sex_info, paste0(sample_name, "_xy_intensities.tsv"), row.names = F, sep = "\t", quote = FALSE)

  
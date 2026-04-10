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
#' sample_info_df          tsv file based on BeadArray sample_sheet with columns Sample_IDAT and Gender
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

sample_info_df <- read.csv("${sample_sheet}", sep = ",", skip = 7)
meth_QC <- readRDS("${meth_rds}")
sample_name <- "${sample_id}"

#### --- Compute X and Y intensities

# Reorder label to ensure they are in the same order in the meth object and in the sample sheet given by the user
sample_info_df[["Sample_IDAT"]] <- factor(sample_info_df[["Sample_IDAT"]], levels = meth_QC[["meta"]][["sample_id"]])
sample_info_df <- sample_info_df[order(sample_info_df[["Sample_IDAT"]]), ]

# Rename column to match predict_sex function output
sample_info_df[["Gender"]][sample_info_df[["Gender"]] == "M"] <- "m"
sample_info_df[["Gender"]][sample_info_df[["Gender"]] == "F"] <- "f"

if (colSums(is.na(meth_QC[["M"]] + meth_QC[["U"]])) <= 1000){
  
  # Apply dye-biaised correction
  meth_extra_QC_corrected <- meth_QC %>% correct_dye_bias 
  # Extract normalized average X and Y intensities
  sex_info <- ewastools::check_sex(meth_extra_QC_corrected)
  
  # Combine it with the sample dataframe 
  sex_info_df <- as.data.frame(sex_info)
  rownames(sex_info_df) <- meth_QC[["meta"]][["sample_id"]]
  sex_info_df[["Sample_IDAT"]] <- meth_QC[["meta"]][["sample_id"]]
  sample_info_df <- sample_info_df %>% left_join(sex_info_df, by = "Sample_IDAT") %>%
    drop_na(Sample_IDAT)
  
  # Keep only relevant information
  sex_info <- sample_info_df %>% select(c(Sample_Name, Sample_IDAT, X, Y, Gender))

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

  
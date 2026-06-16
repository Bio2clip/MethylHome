#!/usr/bin/env Rscript

# Compute QC metrics genome studio like
# Libraries : ewastools, dplyr, tidyr

#### --- DOCUMENTATION --- ####
#' Extract raw values of quality control probes. 
#'
#' @Description 
#' Extract raw values of quality control probes. 
#' It gives intensities for all probe used in genome studio (Illumina).
#'
#' It generates a TSV file with sample_name, sample_idat, channel, metrics and values.
#'
#' @Usage 
#' 
#' compute_qc_gs(meth_rds, sample_id)
#' 
#' @Arguments 
#' 
#' meth_rds      R data file as returned by load_idats     
#' sample_id       Sample name
#'
#### --------------------------- ####

library(ewastools)
library(dplyr)
library(tidyr)

meth_qc <- readRDS("${meth_rds}")
sample_name <- "${sample_id}"

# Store control values in red and green channels
ctrlG <- meth_qc[["ctrlG"]]
ctrlR <- meth_qc[["ctrlR"]]

# Plot Restoration (only green channel)

# Retrieve probes indexes
print(unique(meth_qc[["controls"]][["group"]]))
print(meth_qc[["controls"]][meth_qc[["controls"]][["group"]] == "RESTORATION", index])
index <- meth_qc[["controls"]][meth_qc[["controls"]][["group"]] == "RESTORATION", index]
control_restoration <- as.data.frame(t(ctrlG[index,]))
colnames(control_restoration) <- "value"
control_restoration[["Metric"]] <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == "RESTORATION"]
rownames(control_restoration) <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == "RESTORATION"]
control_restoration[["Global_metrics"]] <- "RESTORATION"
control_restoration[["Sample_Name"]] <- sample_name
control_restoration[["Sample_IDAT"]] <- meth_qc[["meta"]][["sample_id"]]
control_restoration[["Channel"]] <- "green"

# Other wet lab QC metrics
ctrls <- unique(meth_qc[["controls"]][["group"]])
to_remove <- c("NORM_T", "NORM_G", "NORM_C", "NORM_A", "NEGATIVE", "RESTORATION")

# Remove unwanted control metrics
ctrls <- ctrls[!(ctrls %in% to_remove)]

# Create empty dataframe
all_control_green <- control_restoration[-1,]
all_control_red <- control_restoration[-1,]

for (control in ctrls) {
  
  # Red channel
  index <- meth_qc[["controls"]][meth_qc[["controls"]][["group"]] == control, index]
  control_red <- as.data.frame(ctrlR[index,])
  colnames(control_red) <- "value"
  control_red[["Metric"]] <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == control]
  rownames(control_red) <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == control]
  control_red[["Global_metrics"]] <- control
  control_red[["Sample_Name"]] <- sample_name
  control_red[["Sample_IDAT"]] <- meth_qc[["meta"]][["sample_id"]]
  control_red[["Channel"]] <- "red"
  
  
  # Green channel
  control_green <- as.data.frame(ctrlG[index,])
  colnames(control_green) <- "value"
  control_green[["Metric"]] <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == control]
  rownames(control_green) <- meth_qc[["controls"]][["name"]][meth_qc[["controls"]][["group"]] == control]
  control_green[["Global_metrics"]] <- control
  control_green[["Sample_Name"]] <- sample_name
  control_green[["Sample_IDAT"]] <- meth_qc[["meta"]][["sample_id"]]
  control_green[["Channel"]] <- "green"
  
  
  # Merge data frames for green and red separately
  all_control_green <- rbind(all_control_green, control_green)
  all_control_red <- rbind(all_control_red, control_red)
  
}

# Order columns 
control_restoration <- control_restoration[,c("Sample_Name", "Sample_IDAT", "Channel","Global_metrics","Metric","value")]
all_control_green <- all_control_green[,c("Sample_Name", "Sample_IDAT", "Channel","Global_metrics","Metric","value")]
all_control_red <- all_control_red[,c("Sample_Name", "Sample_IDAT", "Channel","Global_metrics","Metric","value")]

# Merge all channel in a single dataframe
all_control <- rbind(control_restoration, all_control_green, all_control_red )

# Write TSV file
write.table(all_control, paste0(sample_name, "_qc_metrics_gs.tsv"), row.names = F, sep = "\t", quote = FALSE, dec = ".")

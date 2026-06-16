#!/usr/bin/env Rscript

# Load IDATs
# Libraries : ewastools, stringr

#### --- DOCUMENTATION --- ####
#' Read and store idat files with minfi package
#' 
#' @Description 
#' 
#'Read idat files and store them in an R object as an RgChannelSet.
#' 
#' @Usage 
#' 
#' load_idats_minfi(files, sample_name)
#' 
#' @Arguments 
#' 
#' files          idat Green files 
#' sample_name    name of the given samples
#' 
#' 

library(minfi)
library(stringr)

### --- Read idat files 

files <- "${idat_green}"
sample_name <- "${sample_id}"

# Get rid of "_Grn.idat" 
infile <- gsub("_Grn.idat", "", x=files)

rgset <- read.metharray(infile) 

### --- Save output for Nextflow
saveRDS(rgset, file = paste0(sample_name, "_minfi.rds"))

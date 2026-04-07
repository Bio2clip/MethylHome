#!/usr/bin/env Rscript

# Load IDATs
# Libraries : ewastools, stringr

#### --- DOCUMENTATION --- ####
#' Read and store idat files 
#' 
#' @Description 
#' 
#'Read idat files and store them in an R object.
#' 
#' @Usage 
#' 
#' load_idats(files, sample_name)
#' 
#' @Arguments 
#' 
#' files          idat Green files 
#' sample_name    name of the given samples
#' 
#' 

library(ewastools)
library(stringr)

### --- Read idat files 

files <- "${idat_green}"
sample_name <- "${sample_id}"

# Get rid of "_Grn.idat" 
infile <- gsub("_Grn.idat", "", x=files)

print(paste("Sample prefix: ", infile))

meth <- read_idats(infile, quiet=FALSE) # `quiet=TRUE` suppresses the progress bar

### --- Save output for Nextflow
saveRDS(meth, file = paste0(sample_name, ".rds"))


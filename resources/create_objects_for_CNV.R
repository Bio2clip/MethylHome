#!/usr/bin/env Rscript


# ----------------------------------------------------------------------
# Script : create_objects_for_CNV.R
# Object  : Create and store objects for reference and annotation to compute CNV
# Librairies : minfi, conumee2
# ----------------------------------------------------------------------


#### --- DOCUMENTATION --- ####
#' Create and store objects for reference and annotation to compute CNV from IDAT files 
#' 
#' @Description 
#' 
#' Create and store objects for reference and annotation to compute CNV from IDAT files from GEO GSE306226.  
#' 
#' @Usage 
#' 
#' Rscript create_objects_for_CNV.R GEO_idats_folder output_directory
#' 
#' @Arguments 
#' 
#' idat_dir      directory of the GSE306226 idat files from the reference database
#' outdir        directory on which to write reference and annotation files
#' 
#' @Value 
#' 
#' Write a .Rdata files for :
#'    * ref_f : reference object for female samples for CNV computation
#'    * ref_m : reference object for male samples for CNV computation
#'    * ref_mf : reference object for unknown sex samples for CNV computation
#'    * anno_epic : annotation object for CNV computation
#'    

library(minfi)
library(conumee2)


# Get filepaths to GEO IDAT files
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("stop Usage : Rscript create_objects_for_CNV.R <GEO_IDAT_directory> <outdir>")
}

GEOdir <- normalizePath(args[1])
outdir <- normalizePath(args[2]) 

# separate male and female files
files_m <- c("GSM9195035_205862290042_R06C01", "GSM9195036_205862290042_R07C01", "GSM9195039_205862290025_R06C01",
             "GSM9195041_205862290025_R08C01", "GSM9195046_205854150036_R05C01", "GSM9195048_205854150036_R07C01",
             "GSM9195052_205854150047_R07C01", "GSM9195053_205854150047_R08C01", "GSM9195054_205854150069_R05C01",
             "GSM9195055_205854150069_R06C01", "GSM9195075_205854150040_R02C01", "GSM9195077_205854150040_R03C01",
             "GSM9195082_205854150049_R01C01", "GSM9195084_205854150049_R03C01", "GSM9195102_205854150145_R01C01",
             "GSM9195104_205854150145_R03C01")
files_f <- c("GSM9195042_205862290011_R05C01", "GSM9195045_205862290011_R08C01", "GSM9195058_205854150020_R05C01",
             "GSM9195059_205854150020_R06C01", "GSM9195062_205854150133_R05C01", "GSM9195063_205854150133_R06C01",
             "GSM9195067_205854150142_R06C01", "GSM9195069_205854150142_R08C01", "GSM9195070_205854150152_R01C01",
             "GSM9195072_205854150152_R03C01", "GSM9195078_205854150051_R05C01", "GSM9195079_205854150051_R06C01",
             "GSM9195086_205854150093_R01C01", "GSM9195089_205854150093_R04C01", "GSM9195091_205854150154_R06C01",
             "GSM9195092_205854150154_R07C01", "GSM9195095_205854150021_R06C01", "GSM9195097_205854150021_R08C01",
             "GSM9195099_205854150146_R06C01", "GSM9195101_205854150146_R08C01", "GSM9195107_206513950148_R02C01",
             "GSM9195108_206513950148_R03C01", "GSM9195112_206513950183_R03C01", "GSM9195113_206513950183_R04C01")

# Create Mset objects for all categories : male, female, mixed
RGsetmodele_m_epic <- read.metharray(paste0(GEOdir, "/", files_m))
RGsetmodele_f_epic <- read.metharray(paste0(GEOdir, "/", files_f))
Msetmodele_m_epic <- preprocessRaw(RGsetmodele_m_epic)
Msetmodele_f_epic <- preprocessRaw(RGsetmodele_f_epic)
Msetmodele_mf_epic <- preprocessRaw(minfi::combine(RGsetmodele_f_epic, RGsetmodele_m_epic))

# Create CNV.load objects for all categories : male, female, mixed
data.cM_epic <- conumee2::CNV.load(Msetmodele_m_epic)
data.cF_epic <- conumee2::CNV.load(Msetmodele_f_epic)
data.cMF_epic <- conumee2::CNV.load(Msetmodele_mf_epic)

# Create annotation for 
anno_epic <- conumee2::CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000,   bin_maxsize = 5000000, array_type = c("EPIC","EPICv2"), exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = TRUE)


save(data.cF_epic, file = file.path(outdir, "data/epic_geo_ref_f.Rdata"))
save(data.cM_epic, file = file.path(outdir, "data/epic_geo_ref_m.Rdata"))
save(data.cMF_epic, file = file.path(outdir, "data/epic_geo_ref_mf.Rdata"))
save(anno_epic, file = file.path(outdir, "data/annoXY_epic.Rdata"))
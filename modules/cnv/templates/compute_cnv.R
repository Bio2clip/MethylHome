#!/usr/bin/env Rscript

# Compute cnv from methylation data (IDATs)
# Libraries : conumee2, minfi, IlluminaHumanMethylationEPICv2manifest, IlluminaHumanMethylationEPICv2anno.20a1.hg38

#### --- DOCUMENTATION --- ####
#' Compute cnv from methylation data and generates plots and tables.
#' 
#' @Description 
#' 
#'Compute cnv from methylation data, write data in tsv, seg and IGV files, generates visualizatuions at different levels (genome, chr, gene).
#'The input file is based on BeadArray sample_sheet with columns Sample_IDAT and Gender
#' 
#' @Usage 
#' 
#' CNV(ref_f, ref_m, anno, sample_id, idat, sample_sheet)
#' 
#' @Arguments 
#' 
#' ref_f                   R data object containing reference dataset for female samples
#' ref_m                   R data object containing reference dataset for male samples
#' anno                    R data object containing annotation for reference in EPICv1 and query in EPICv2 and details set to DKFZ list
#' anno_path_cosmic        R data object containing annotation for reference in EPICv1 and query in EPICv2 and details set to cosmic genes
#' sample_id               name of the given samples
#' idat                    IDAt file of the query sample
#' sex                     sex of the studied sample
#' CNV_focal               TRUE or FALSE, wether to compute CNV.focal()
#' 
#' @Note
#' Because of CNV.focal(), CNV.detailplot() and CNV.detailplot_wrap() functions, an internet connection is required.
#' An option to disable these execution is implemented : CNV.focal = FALSE

#### --------------------------- ####


library(conumee2)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(stringr)

# Retrieve reference and annotation objects
refF_path <- "$ref_f"
refM_path <- "$ref_m"
refMF_path <- "$ref_mf"

anno_path <- "$anno"
anno_path_cosmic <- "$anno_cosmic"
CNV_focal <- "$CNV_focal"

# Get query sample name and sex
sample_name <- "${sample_id}"
sex <- "${sex}"


# Load data ; rgsetobject
RGsettarget <- readRDS("${rgset_rds}")
sample_IDAT <- RGsettarget@colData@rownames[1]

# Dupplicate file with dumy name to have more than one raw in the rgset object (to make conumee2 work)
RGsetAdd <- RGsettarget
RGsetAdd@colData@rownames <- "dumy"

RGset <- minfi::combine(RGsettarget, RGsetAdd)
Msettarget <- minfi::preprocessRaw(RGset)

# Select reference dataset depending on sample sex
if(sex == "F" || sex == "Female"){
  RGset <- minfi::combine(RGsettarget, RGsetAdd)
  Msettarget <- minfi::preprocessRaw(RGset)
  load(refF_path)
  ref <- data.cF_epic
} else if (sex == "M" || sex == "Male") {

  load(refM_path)
  ref <- data.cM_epic
  
} else {

    load(refMF_path)
    ref <- data.cMF_epic
}

# Load annotation object
load(anno_path)

# Create cnv query object
query <- conumee2::CNV.load(Msettarget)

# Compute CNV
x <- conumee2::CNV.fit(query = query, ref = ref, anno_epic)
x <- conumee2::CNV.bin(x)
x <- conumee2::CNV.detail(x)
x <- conumee2::CNV.segment(x)

if (CNV_focal == "true") {
  x <- conumee2::CNV.focal(x, sig_cgenes = F)
}


# Retrieve sample idat
SampleIDAT<-RGsettarget@colData@rownames[1]

# Global CNV plot
png(paste0(sample_name, "_", SampleIDAT,"_AllChr.png"), width = 40, height = 20, units = 'cm', res = 300)
conumee2::CNV.genomeplot(x[1])
dev.off()

# CNV plot by chr
for (i in 1:22) {
  chr<-paste0("chr",i)
  png(paste0(sample_name, "_", SampleIDAT,"_",chr,".png"), width = 20, height = 15, units = 'cm', res = 300)
  conumee2::CNV.genomeplot(x[1], chr = chr)
  dev.off()
}
print(CNV_focal)
# CNV plot by detail regions
if (CNV_focal == "true") {
  genes<- anno_epic@detail\$name # Retrieve detail regions names
  print(genes)
  for (i in 1:length(genes)) {
    gene<-genes[i]
    print(gene)
    png(paste0(sample_name, "_", SampleIDAT, "_",gsub(x=gene, "/", "-"),"_gene.png"), width = 20, height = 15, units = 'cm', res = 300)
    conumee2::CNV.detailplot(x[1], name = gene)
    dev.off()
}

# Global detail regions CNV plot
png(paste0(sample_name, "_", SampleIDAT, "_Genes",".png"), width = 40, height = 20, units = 'cm', res = 300)
conumee2::CNV.detailplot_wrap(x[1])
dev.off()
}


# Compute general metrics from cnv computation
metrics <- data.frame(Sample_Name = sample_name,
                      Sample_IDAT = SampleIDAT,
                      noise = x[1]@fit[["noise"]],
                      shift = x[1]@bin[["shift"]])


for (i in rownames(x[1]@fit[["coef"]]) ) {
  metrics[[i]] <- x[1]@fit[["coef"]][i,]
}

# Write igv and text files 
segments <- conumee2::CNV.write(x[1], what = "segments") 
segments[[1]][["Sample_Name"]] <- sample_name
write.table(segments[[1]], paste0(sample_name, "_", SampleIDAT,"_CNVsegments.seg"), sep = "\t", quote = F, row.names = F)

detail <- conumee2::CNV.write(x[1], what = "detail")
detail[[1]][["Sample_Name"]] <- sample_name
write.table(detail[[1]], paste0(sample_name, "_", SampleIDAT, "_CNVdetail.tsv"), sep = "\t", quote = F, row.names = F)

conumee2::CNV.write(x[1], what = "bins", file = paste0(sample_name, "_", SampleIDAT, "_CNVbins.igv"))
conumee2::CNV.write(x[1], what = "probes", file = paste0(sample_name, "_", SampleIDAT, "_CNVprobes.igv"))

write.table(metrics, paste0(sample_name, "_", SampleIDAT, "_metrics.tsv"), sep = "\t", row.names = F, quote = F)


### Compute values for all cosmic genes

load(anno_path_cosmic)

# Compute CNV
y <- conumee2::CNV.fit(query = query, ref = ref, anno_epic_cosmic)
y <- conumee2::CNV.bin(y)
y <- conumee2::CNV.detail(y)
y <- conumee2::CNV.segment(y)

if (CNV_focal == "true") {
  y <- conumee2::CNV.focal(y, sig_cgenes = F)
}

# Retrieve sample idat
SampleIDAT<-RGsettarget@colData@rownames[1]


detail_cosmic <- conumee2::CNV.write(y[1], what = "detail")
detail_cosmic[[1]][["Sample_Name"]] <- sample_name
write.table(detail_cosmic[[1]], paste0(sample_name, "_", SampleIDAT, "_CNVdetail_cosmic.tsv"), sep = "\t", quote = F, row.names = F)

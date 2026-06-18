#!/usr/bin/env Rscript

#### --- DOCUMENTATION --- ####
#' Compute tumor purity of the given sample
#' 
#' @Description 
#' 
#' Computes tumor purity of the given sample. No normalization is performed.
#' 
#' @Usage 
#' 
#' compute_tumor_purity(rgset_rds, sample_id)
#' 
#' @Arguments 
#' 
#' rgset_rds      Output of calling minfi::read_metharray.
#' sample_id      Name of the sample
#' 
#' @Value 
#' 
#' Return a dataframe containing all computed metrics per sample.
#' 
#' @Note 
#' 
#' Since RFpurify consists of RAndom forest models trained on 450k Beadchhip, the script converts EPICv2 probes names to EPIC names. 
#' Some probes are missing on EPICv2 Beadchip: 0.5 default value is put for these probes. 

# Load packages 
library(RFpurify)
library(minfi)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)

# Get query sample name and file
sample_name <- "${sample_id}"
query_rgset <- readRDS("${rgset_rds}")

Mset_m <- preprocessRaw(query_rgset)

if (sum(getUnmeth(Mset_m) == 0 | getMeth(Mset_m) == 0) < 1000) {
  # Compute beta values
  betas_epicv2 <- getBeta(Mset_m)

  # Transform probe names to fit previous annotation
  rownames(betas_epicv2) <-  gsub("_.*\$", "", x=rownames(betas_epicv2))

  # Load model to retrieve probes of interest
  data("RFpurify_ABSOLUTE")  
  data("RFpurify_ESTIMATE") 
  model_probes_absolute <- rownames(RFpurify_ABSOLUTE\$importance)
  model_probes_estimates <- rownames(RFpurify_ESTIMATE\$importance)

  # Retrieve probes only present in 450k
  diff_absolute <- model_probes_absolute[!model_probes_absolute %in% rownames(betas_epicv2) ]
  diff_estimates <- model_probes_estimates[!model_probes_estimates %in% rownames(betas_epicv2) ]

  # Missing probes matrix with rownames to probes id
  mat_test_abs <- matrix(
    rep(0.5, ncol(betas_epicv2)),
    nrow = length(diff_absolute),
    ncol = ncol(betas_epicv2),
    byrow = FALSE
  )
    
  rownames(mat_test_abs) <- diff_absolute
    
  mat_test_est <- matrix(
    rep(0.5, ncol(betas_epicv2)),
    nrow = length(diff_estimates),
    ncol = ncol(betas_epicv2),
    byrow = FALSE
  )
  rownames(mat_test_est) <- diff_estimates
    
  # Complete beta values matrix
  betas_allcpg_v2 <- rbind(betas_epicv2, mat_test_abs, mat_test_est)
    
  # Keep only probes of interest
  betas_clean_cpg_v2 <- betas_allcpg_v2[rownames(betas_allcpg_v2) %in% c(model_probes_estimates, model_probes_absolute),, drop=F]
    
  # Keep unique values for duplicated probes
  dup <- rownames(betas_clean_cpg_v2)[duplicated(rownames(betas_clean_cpg_v2))]
    
  for (cpg in dup){
    betas_clean_cpg_v2[rownames(betas_clean_cpg_v2) == cpg]  <- mean(betas_clean_cpg_v2[rownames(betas_clean_cpg_v2) == cpg] )
  }
    
  betas_m <- betas_clean_cpg_v2[!duplicated(rownames(betas_clean_cpg_v2)),, drop=FALSE]
    
  # Predict tumor purity with ABSOLUTE and ESTIMATE method
  absolute <- predict_purity_betas(betas_m,method="ABSOLUTE")
  estimate <- predict_purity_betas(betas_m,method="ESTIMATE")

  # Create dataframe with computed purities
  purity <- data.frame(Sample_Name = Mset_m@colData@rownames,
                        absolute = absolute,
                        estimate = estimate)
    
  # Create table to display
  tumor_purity_df <- data.frame(analysis = c("estimate", "absolute"),
                    value = c(purity\$estimate, purity\$absolute)
                    )

  custom_theme <- ttheme_default(
    core = list(base_size = 16))

  final_purity_table <- tableGrob(tumor_purity_df, rows = NULL, theme = custom_theme)

  # Estimate - Absolute plots
  purity_plot <- ggplot(purity, aes(x = estimate, y = absolute)) + 
    geom_point(size = 4) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray36") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray36") + 
    theme_minimal()
} else {
  purity_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Too low quality to display \n tumor purity", size = 5, hjust = 0.5, col="red") +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
  
  final_purity_table <- NULL
  purity <- data.frame(Sample_Name = Mset_m@colData@rownames,
                        absolute = NA,
                        estimate = NA)
}

# Write PDF file
pdf(paste0(sample_name, "_tumor_purity.pdf"),  width = 14, height = 10)
(purity_plot| final_purity_table) + plot_annotation(title = paste0("Tumor purity Report - ", sample_name ),
                                      subtitle = paste("Date :", Sys.Date(), "; Library: RFpurify"),
                                      theme = theme(plot.title = element_text(size = 20, face = "bold"))) +
  plot_layout(widths = c(1.4, 0.6))
dev.off()

write.table(purity, paste0(sample_name,"_tumor_purity.tsv"), sep = "\t", row.names = F, quote = F)

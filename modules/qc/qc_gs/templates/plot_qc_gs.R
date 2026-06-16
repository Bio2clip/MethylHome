#!/usr/bin/env Rscript

# Plot QC metrics genome studio like
# Libraries : ewastools, dplyr, tidyr, ggplot2, patchwork

#### --- DOCUMENTATION --- ####
#' Plot raw values of quality control probes. 
#'
#' @Description 
#' Extract raw values of quality control probes. 
#' It generates plots for all metrics like in genome studio (Illumina).
#'
#'
#' @Usage 
#' 
#' plot_qc_gs(control_metrics_gs)
#' 
#' @Arguments 
#' 
#' control_metrics_gs      TSV file containing values for all metrics values (including restoration) in both red and green chanel.    
#'

library(ewastools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

### --- Plots

control_all <- read.csv("$control_metrics_gs", sep = "\t")

# Split channels by color and metrics
control_red <- control_all %>% filter(Channel == "red")
control_green <- control_all %>% filter(Channel == "green", Global_metrics != "RESTORATION")
control_restoration <- control_all %>% filter(Global_metrics == "RESTORATION")

# Create color pallet for metrics plots
pallet <- c("red", "dodgerblue2", "green4", "darkorchid3", "orange", "yellow", "chocolate4", "pink2","grey","cadetblue2", "coral2", "darkolivegreen2" )

# Restoration
# Format long for ggplot
control_restoration_long <- control_restoration %>%
  pivot_longer(cols = -c( Metric, Global_metrics, Sample_Name, Sample_IDAT, Channel),
               values_to = "Intensity" )

# Plot for restoration control (only in green channel)
restoration_plot <- ggplot(control_restoration_long, aes(x = Sample_Name, y = Intensity, color = Metric)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "RESTORATION GREEN Efficiency",
       x = "Sample",
       y = "Intensity") + 
  scale_color_manual(values = "green4")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("genome_studio_like_plot.pdf", width = 14, height = 10)
print(restoration_plot)

ctrls <- unique(control_red[["Global_metrics"]])

for (control in ctrls) {
  
  # Green channel metrics plot
  control_green2 <- control_green %>% filter(Global_metrics == control)
  
  control_green_long <- control_green2 %>%
    tibble::rownames_to_column("Control") %>%
    pivot_longer(cols = -c(Control, Metric, Global_metrics, Sample_Name, Sample_IDAT, Channel),
                 names_to = "Sample",
                 values_to = "Intensity" )
  
  control_green_plot <- ggplot(control_green_long, aes(x = Sample_Name, y = Intensity, color = Metric)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = paste0(control," GREEN Efficiency"),
         x = "Sample",
         y = "Intensity") + 
    scale_color_manual(values = pallet) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Red channel metrics plot
  control_red2 <- control_red %>% filter(Global_metrics == control)
  
  control_red_long <- control_red2 %>%
    pivot_longer(cols = -c(Metric, Global_metrics, Sample_Name, Sample_IDAT, Channel),
                 values_to = "Intensity" )
  
  control_red_plot <- ggplot(control_red_long, aes(x = Sample_Name, y = Intensity, color = Metric)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = paste0(control, " RED Efficiency"),
         x = "Sample",
         y = "Intensity") + 
    scale_color_manual(values = pallet) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  # Plot
  print((control_green_plot | control_red_plot))
  
}

dev.off()
#!/usr/bin/env Rscript

# Plot QC 
# Libraries : ggplot2, dplyr, ggridges, patchwork, gridExtra, minfi (beta)

#### --- DOCUMENTATION --- ####
#' Display Quality Control 
#' 
#' @Description 
#' 
#' Reproduces plots as displayed in Epignostix (DKFZ). It produces a plot that consists of a
#' QC report. Plots displayed are : 
#'  - 21 metrics described in BeadArray Software Guide (Illumina),
#'  - log2 methylated and unmethylated intensities,
#'  - ECDF detection p-value (of the selected sample),
#'  - beta value distribution per probe type (of the selected sample).
#' 
#' Plots are saved in a pdf file.
#' 
#' @Usage 
#' 
#' plot_qc(qc_tsv, database, meth_rds, sample_id)
#' 
#' @Arguments 
#' 
#' qc_tsv          tsv file with qc metrics as returned by extract_qc_metrics
#' database        tsv file consisting of the output of extract_qc_metrics for the reference idat files
#' meth_rds        R data file as returned by load_idats
#' sample_id       Sample name
#' 
#' 
 
#### --- LIBRAIRIES --- ####

library(dplyr)
library(tidyr)
library(ewastools)
library(ggplot2)
library(ggridges)
library(patchwork)
library(gridExtra)

### --- INPUT --- ###

### --- Retrieve qc metrics
qc_df <- read.csv("${qc_tsv}", sep = "\t")
colnames(qc_df) <- gsub("_", " ", x=colnames(qc_df))

### --- Retrieve database
database <- read.csv("${database}", sep = "\t")
colnames(database) <- gsub("_", " ", x=colnames(database))

### --- Retrieve idats object
meth_QC <- readRDS("${meth_rds}")

sample_id <- qc_df[["Sample"]]
rownames(qc_df) <- qc_df[["Sample"]]
sample_name <- "${sample_id}"
  
### --- Log2Meth and unmeth intensities quality control
  
log2_un_meth_intensities_plot <- ggplot(database, aes(x = Log2UnmethIntensity, y = Log2MethIntensity)) +
    geom_point(alpha = 0.7) + 
    geom_vline(xintercept = 8, linetype = "dashed", color = "red") + # vertical line to x=8
    geom_hline(yintercept = 8, linetype = "dashed", color = "red") + # horizontal line to y=8
    geom_point(data = qc_df, 
               aes(x = as.numeric(Log2UnmethIntensity), y = as.numeric(Log2MethIntensity)), 
               color = "red", size = 2) +
    coord_cartesian(xlim = c(7.5, 12.5), ylim = c(7.5, 12.5)) +
    theme_minimal() + 
    labs(
      x = "Median Log2 Unmethylated Intensity",
      y = "Median Log2 Methylated Intensity",
      title = "QC without Normalization"
    ) + theme(plot.title = element_text(size=11, face="bold", vjust = 0.5, hjust = 0.5))

### --- Plot ecdf detP
detP_ewas <- meth_QC %>% detectionP
detP_df <- as.data.frame(detP_ewas[["detP"]]) # Convert to df for ggplot
colnames(detP_df) <- meth_QC[["meta"]][["sample_id"]]
ecdf_P <- ecdf(detP_df[[sample_id]]) # Compute ecdf of the studied sample

  # Plot
ecdf_detP_plot <- ggplot(detP_df, aes(.data[[sample_id]])) + stat_ecdf(geom = "step")+
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = 0.05, 
      y = 0.05,         
      label = paste0("DetectionP < 0.05 : ", round(ecdf_P(0.05), 3)),
      color = "red",
      hjust = -0.1
    ) +
    labs(title = "Detection Rate", y = "ECDF", x="Detection p-values")+
    theme_minimal() + 
    theme(plot.title = element_text(size=11, face="bold", vjust = 0.5, hjust = 0.5))
  
### --- Quality metrics BeadArray
table_res_t <- database %>% select(-c(Sample, DetectionRate, Log2MethIntensity, Log2UnmethIntensity))

# Reference dataset for density distribution
df_long <- table_res_t %>%
    mutate(SampleID = rownames(table_res_t)) %>%
    pivot_longer(
      cols = -SampleID,
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(metric) %>% # group by to compute sd / metric and not a global sd
    mutate(
      sd_metric = sd(value, na.rm = TRUE),
      value = value / sd_metric # Reduction by QC metrics sd
    ) %>%
    ungroup()

# Compute thresholds and scaled thresholds

df_threshold <- data.frame(
    metric = colnames(table_res_t),
    threshold = as.numeric(c(0, 5,5,5,5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5))
  ) %>%
    left_join(
      df_long %>% group_by(metric)) %>% select(-c(SampleID, value))
df_threshold[["threshold_scaled"]] <- df_threshold[["threshold"]] / df_threshold[["sd_metric"]]
df_threshold <- unique(df_threshold)

# Create Sample dataframe
df_sample <- as.data.frame(t(qc_df %>% select(-c(DetectionRate, Log2MethIntensity, Log2UnmethIntensity)))) 
colnames(df_sample) <- as.character(df_sample[2,])
df_sample <-  df_sample[-c(1,2), , drop = FALSE] %>% select(all_of(sample_id))

df_sample <- df_sample %>%
  dplyr::select(all_of(sample_id))
colnames(df_sample) <- "value"
df_sample[["SampleID"]] <- sample_id
df_sample[["metric"]] <- rownames(df_sample)
df_sample <- df_sample  %>%
  left_join(
    df_threshold, by = "metric") %>% select(-c(threshold, threshold_scaled))
df_sample[["scaled_value"]] <- as.numeric(df_sample[["value"]]) / as.numeric(df_sample[["sd_metric"]])

# Order label
df_long[["metric"]] <- factor(df_long[["metric"]],
                         levels = rev(colnames(table_res_t)))
df_threshold[["metric"]] <- factor(df_threshold[["metric"]],
                              levels = levels(df_long[["metric"]]))

df_threshold[["y"]] <- as.numeric(df_threshold[["metric"]])

# Create plot
p <-ggplot(df_long, aes(x = value, y = metric)) +

  geom_density_ridges(rel_min_height = 0.005, stat = "density_ridges", fill = "grey85", color = "black", scale = 1) + #rel_min_height argument to remove tails
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_segment(
    data = df_threshold,
    aes(x = threshold_scaled,
        xend = threshold_scaled,
        y = y - 0.4,
        yend = y + 0.4),
    color = "red",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +

  geom_point(
    data = df_sample,
    aes(x = scaled_value, y = metric),
    color = "red",
    size = 2
  ) +

  theme_minimal() +
  labs(x = "Values scaled by SD", y = "Array quality metrics") +
  theme(axis.text.x = element_blank(),  # Get rid of scale x numbers
        axis.ticks.x = element_blank())  # Get rid of scale x ticks
  
  
### --- Plot beta
# Compute beta 
beta = as.data.frame(meth_QC %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize)

# Create dataframe manifest info
reduced_manifest <- meth_QC[["manifest"]] %>%
  select(ilmn_id, probe_design)
  
# To join df
beta[["ilmn_id"]] <- rownames(beta)
  
# Add probe type annotations
beta <- beta %>%
  left_join(reduced_manifest, by = "ilmn_id") %>%
  mutate(Type = probe_design) %>%
  select(-probe_design)
  
# Prepare data for ggplot2
beta_long <- beta %>%
  select(all_of(sample_id), Type) %>%
  filter(Type %in% c("I", "II")) %>%
  rename(Beta = all_of(sample_id)) %>%
  drop_na() %>%
  mutate(Type = factor(Type, levels = c("I", "II")))

# Plot
beta_plot <- ggplot(beta_long, aes(x = Beta, fill = Type, color = Type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 alpha = 0.4,
                 position = "identity") +
  # geom_density(alpha = 0.8, linewidth = 1) +
  scale_fill_manual(values = c("I" = "gray41",
                               "II" = "black")) +
  scale_color_manual(values = c("I" = "gray41",
                                "II" = "black")) +
  labs(title = paste("Distribution per Probe Type"),
       x = "Beta",
       y = "Density",
       fill = "Probe Type",
       color = "Probe Type") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=11, face="bold", vjust = 0.5, hjust = 0.5))

### --- Create summary table / sample

qc_sample_df <- qc_df %>%
  filter(Sample %in% sample_id)

qc_sample_df <- as.data.frame(t(qc_sample_df))
colnames(qc_sample_df) <- "Value"
qc_sample_df <- qc_sample_df[-c(1,2), , drop=F]
qc_sample_df[["Cutoff"]] <- as.numeric(c(0.95, 8, 8, 0, 5,5,5,5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5))
qc_sample_df[["Status"]] <- ifelse(as.numeric(qc_sample_df[["Value"]]) > as.numeric(qc_sample_df[["Cutoff"]]), "pass", "fail")
qc_sample_df[["Value"]] <- round(as.numeric(qc_sample_df[["Value"]]), 3)

# Color according pass/fail status
row_fill <- case_when(
  qc_sample_df[["Status"]] == "fail" ~ "red2",
  qc_sample_df[["Status"]] == "pass"  ~ "gray95")

# Create the theme
custom_theme <- ttheme_default(
  core = list(
    bg_params = list(fill = row_fill)  ),
  base_size = 9
)

table_plot <- tableGrob(qc_sample_df, theme = custom_theme) # Convert to a graphic object for it to be ploted
  
### --- Create plot/ QC Report
final_plot <- (p|(log2_un_meth_intensities_plot/ecdf_detP_plot/beta_plot)|table_plot) + plot_layout(widths = c(1, 0.9, 1.1)) + 
  plot_annotation(
    title = paste0("QC Report - Sample ", sample_name),
    subtitle = paste0("Date:", Sys.Date(), "; Library: ewastools"),
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )
  
### --- Generate pdf report
output_file <- paste0(sample_name, "_qc_plot.pdf")
pdf(output_file, width = 14, height = 8)
print(final_plot)
dev.off()
  

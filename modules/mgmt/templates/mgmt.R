#!/usr/bin/env Rscript

#### --- DOCUMENTATION --- ####
#' Compute MGMT methylation status 
#' 
#' @Description 
#' 
#' Computes MGMT methylation status based on probes cg12434587 and cg12981137. No normalization is applied. 
#' 
#' @Usage 
#' 
#' mgmt(rgset_rds, sample_id)
#' 
#' @Arguments 
#' 
#' rgset_rds      Output of calling minfi::read_metharray.
#' sample_id      Name of the sample
#' 
#' @Value 
#' 
#' Return a pdf file containing a plot for mgmt methylation value and a table with confidence intervals.

require(mgmtstp27)
library(minfi)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

# Get query sample name and file
sample_name <- "${sample_id}"
query_rgset <- readRDS("${rgset_rds}")

dat <- preprocessRaw(query_rgset)
  
sample_id <- query_rgset@colData@rownames

# Compute M-value
mvalue <- as.data.frame(t(getM(dat)))
  
# Predict MGMT status
pred1 <- MGMTpredict(mvalue, ic.distrib = "normal")
  
# quality control graphics
# par(mfrow=c(2,3))
# MGMTqc.pop(pred1,which.plot=1:3,mfrow=NULL)
# MGMTqc.single(pred1,nsample=1,which.plot=1:3,mfrow=NULL)
  
pred1\$state[pred1\$state == "M"] <- "Methylated"
pred1\$state[pred1\$state == "U"] <- "Unmethylated"

# filter to obtain only necessary info 
prediction <- pred1 %>%
  dplyr::select(state, pred, lower, upper, cg12434587, cg12981137)

betas <- getBeta(dat)

prediction[["cg12434587_beta-value"]] <- betas["cg12434587_BC11",]
prediction[["cg12981137_beta-value"]] <- betas["cg12981137_TC11",]

prediction <- prediction %>% 
    mutate(across(where(is.numeric), ~ round(., 3)))
    
custom_theme <- ttheme_default(
  core = list(base_size = 14))

# For the table to be displayed
table_mgmt_plot <- tableGrob(t(prediction), theme = custom_theme)
  
# Plot prediction
mgmt_plot <- ggplot(prediction, aes(y = 0, x = pred, xmin = lower, xmax = upper)) +
  geom_pointrange(color = "black", size = 1) +
  xlim(0,1) +
  geom_vline(xintercept = 0.3582, linetype = "dashed", color = "red") + 
  labs(x = "Predicted Value") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 

# Plot with figure and table
final_mgmt <- (mgmt_plot | table_mgmt_plot) + 
  plot_annotation(
    title = paste0("MGMT Promotor status prediction - Echantillon ", sample_name),
    subtitle = paste("Date :", Sys.Date(), "; Library: mgmtstp27"),
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  ) +
  plot_layout(widths = c(1.4, 0.6))
  
pdf(paste0("MGMT_plot_minfi_", sample_name, ".pdf"), , width = 14, height = 4)
print(final_mgmt)
dev.off()

prediction["Sample_Name"] <- sample_name
write.table(prediction, paste0(sample_name,"_mgmt.tsv"), sep = "\t", row.names = F, quote = F)

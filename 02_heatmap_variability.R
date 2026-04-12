# DNA Methylation Heatmap Pipeline
# =============================================
# Project: B cell lymphomas DNA methylation analysis
# April 2026

library(ComplexHeatmap)

# Load your data first
load("your_betas.RData")

# Compute SD for rows
betas_file$sd <- apply(betas_file, 1, function(x) sd(x))

# Sort CpGs by decreasing variability (most variable first)
betas_file <- betas_file[order(-betas_file$sd), ]

# Create PDF output
pdf("address/Heatmap.pdf", width = 14)

# Generate heatmap using top 15000 variable CpGs
print(Heatmap(as.matrix(betas_file[1:15000, 1:(ncol(betas_file)-1)]), 
              border = FALSE,                    
              show_row_dend = FALSE,             
              show_column_dend = TRUE,           
              row_title = NULL, 
              column_title = NULL,
              cluster_row_slices = FALSE,        
              cluster_column_slices = FALSE,     
              show_row_names = FALSE,            
              show_column_names = FALSE,         
              show_heatmap_legend = TRUE,        
              cluster_rows = TRUE))             

dev.off()

# DMRcate Pipeline: CLL vs MCL Differential Methylation
# =============================================
# Project: B cell lymphomas DNA methylation analysis
# April 2026

library(dplyr)
library(tidyverse)
library(minfi)
library(DMRcate)

# load("your_data.RData")  # must contain: betas_final, CLL_MCL (with Sample_ID and diagnosis columns)

# -----------------------------------------------------------------------------
# 1. Prepare Beta Matrix
# -----------------------------------------------------------------------------

# Subset to samples present in metadata
betas_clean <- betas_final %>%
  select(all_of(intersect(colnames(betas_final), CLL_MCL$Sample_ID)))

# Convert to M-values and remove SNP-associated probes
M_values      <- minfi::logit2(betas_clean)
Mvalues_noSNP <- rmSNPandCH(M_values, rmcrosshyb = FALSE)

# -----------------------------------------------------------------------------
# 2. Define Groups and Design Matrix
# -----------------------------------------------------------------------------

# Rename columns to diagnosis labels
colnames(Mvalues_noSNP) <- CLL_MCL$diagnosis[match(colnames(Mvalues_noSNP), CLL_MCL$Sample_ID)]

# Extract group factor (CLL as reference, MCL as comparison)
type <- gsub("_.*", "", colnames(Mvalues_noSNP))
type <- factor(type, levels=c("CLL", "MCL"))
design <- model.matrix(~type)   # coef=2 tests MCL - CLL

# -----------------------------------------------------------------------------
# 3. CpG Annotation and DMR Identification
# -----------------------------------------------------------------------------

# Annotate CpGs with differential methylation statistics
myannotation <- cpg.annotate("array", object = Mvalues_noSNP, what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design, coef = 2, fdr = 0.001)

# Identify DMRs using kernel smoothing
dmrcoutput     <- dmrcate(myannotation, lambda = 1000, C = 2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

# -----------------------------------------------------------------------------
# 4. Export Results
# -----------------------------------------------------------------------------

result.df <- as.data.frame(results.ranges)

# meandiff = mean(MCL) - mean(CLL)
# meandiff > 0: hypermethylated in MCL | meandiff < 0: hypomethylated in MCL
cat("Hypermethylated in MCL:", sum(result.df$meandiff > 0), "\n")
cat("Hypomethylated in MCL:", sum(result.df$meandiff < 0), "\n")

write_csv(result.df, "CLL_MCL_DMRs.csv")

# -----------------------------------------------------------------------------
# 5. DMR Visualization
# -----------------------------------------------------------------------------

groups <- c(CLL = "forestgreen", MCL = "orange")
cols <- unname(groups[type])

png("DMR_CLL_MCL.png", width = 1200, height = 800, res = 150)
DMR.plot(ranges = results.ranges, dmr = 3, CpGs = myannotation,
         what = "Beta", arraytype = "EPICv2", phen.col = cols, genome = "hg38")
dev.off()

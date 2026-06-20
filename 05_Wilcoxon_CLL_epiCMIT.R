# epiCMIT Pairwise Wilcoxon Test: CLL Subtypes (m-CLL vs i-CLL vs n-CLL)
# =============================================
# Project: B cell lymphomas DNA methylation analysis
# June 2026

library(dplyr)
library(ggplot2)
library(ggsignif)

# load("your_data.RData")  # must contain: df (with diagnosis and epiCMIT columns)
# -----------------------------------------------------------------------------
# 1. Prepare Diagnosis Factor
# -----------------------------------------------------------------------------
# Set factor levels BEFORE running the test, so the output matrix's row/col
# names match the order we expect when indexing it below
df$diagnosis <- factor(
  df$diagnosis,
  levels = c("m-CLL", "i-CLL", "n-CLL")
)
# -----------------------------------------------------------------------------
# 2. Pairwise Wilcoxon Test
# -----------------------------------------------------------------------------
# Holm correction (default) for multiple comparisons across the three groups
pw <- pairwise.wilcox.test(
  df$epiCMIT,
  df$diagnosis,
  p.adjust.method = "holm"
)
print(pw)
# pairwise.wilcox.test only returns the lower triangle: each pair is indexed
# [later_level, earlier_level] given the factor level order set in step 1
pvec <- c(
  pw$p.value["i-CLL", "m-CLL"],
  pw$p.value["n-CLL", "i-CLL"],
  pw$p.value["n-CLL", "m-CLL"]
)
p_labels <- formatC(pvec, digits = 2, format = "g")
cat("Adjusted p-values:\n")
print(pvec)
cat("Formatted labels:\n")
print(p_labels)
# -----------------------------------------------------------------------------
# 3. Export Results
# -----------------------------------------------------------------------------
result_df <- data.frame(
  comparison = c("m-CLL vs i-CLL", "n-CLL vs i-CLL", "n-CLL vs m-CLL"),
  p_adjusted = pvec,
  p_label    = p_labels
)
write.csv(result_df, "CLL_epiCMIT_wilcoxon_pvalues.csv", row.names = FALSE)
# -----------------------------------------------------------------------------
# 4. Boxplot Visualization with Significance Brackets
# -----------------------------------------------------------------------------
all_pairs <- list(
  c("m-CLL", "i-CLL"),
  c("n-CLL", "i-CLL"),
  c("n-CLL", "m-CLL")
)
y_max <- max(df$epiCMIT, na.rm = TRUE)
p <- ggplot(df, aes(x = diagnosis, y = epiCMIT)) +
  geom_boxplot(fill = "#393b79", alpha = 0.7) +
  geom_signif(
    comparisons = all_pairs,
    annotations = p_labels,
    textsize = 4,
    y_position = c(1.05, 1.10, 1.15) * y_max,
    tip_length = 0.02
  ) +
  labs(
    x = "Diagnosis", y = "epiCMIT",
    title = "epiCMIT by CLL subtype (Wilcoxon)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CLL_wilcoxon_epiCMIT.pdf", plot = p, width = 12, height = 5)

# Wilcoxon pairwise comparisons of epiCMIT across CLL subtypes
library(ggplot2)
library(ggsignif)

# 1. Set factor levels BEFORE running the test, so output matrix dimnames match what we index below
df$diagnosis <- factor(
  df$diagnosis,
  levels = c("m-CLL", "i-CLL", "n-CLL")
)

# 2. Pairwise Wilcoxon test
pw <- pairwise.wilcox.test(
  df$epiCMIT,
  df$diagnosis,
  p.adjust.method = "holm"
)
print(pw)

# 3. Extract the three pairwise p-values
pvec <- c(
  pw$p.value["i-CLL", "m-CLL"],
  pw$p.value["n-CLL", "i-CLL"],
  pw$p.value["n-CLL", "m-CLL"]
)
p_labels <- formatC(pvec, digits = 2, format = "g")
print(pvec)
print(p_labels)

# 4. Plot
all_pairs <- list(
  c("m-CLL", "i-CLL"),
  c("n-CLL", "i-CLL"),
  c("n-CLL", "m-CLL")
)

y_max <- max(CLL_epiCMIT$epiCMIT, na.rm = TRUE)

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

# 5. Save outputs
ggsave("CLL_Wilcoxon_epiCMIT.pdf", plot = p, width = 12, height = 5)

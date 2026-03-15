getwd()

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(knitr) # markdown tables

# Data
DATA_DIR <- "./data/output"
VIZ_DIR <- "./output"

data <- read_tsv(paste0(DATA_DIR, "/data1.tsv"), col_names = TRUE)
head(data)

type_colours <- c(SNP = "skyblue2", INDEL = "olivedrab4")


# Preprocessing
dp_summary <- data %>%
  summarise(
    mean_DP = mean(DP, na.rm = TRUE),
    max_threshold_99th = quantile(DP, 0.99, na.rm = TRUE),
    min_threshold_1st = quantile(DP, 0.01, na.rm = TRUE),
    mean_QUAL = mean(QUAL, na.rm = TRUE),
    max_threshold_99th_QUAL = quantile(QUAL, 0.99, na.rm = TRUE),
    min_threshold_1st_QUAL = quantile(QUAL, 0.01, na.rm = TRUE)
  )

print(dp_summary)


p1 <- ggplot(data, aes(x = DP)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(dp_summary$min_threshold_1st, dp_summary$max_threshold_99th), 
             color = "red", linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Read Depth - raw data",
    x     = "DP (log10 scale)",
    y     = "Count"
  ) + theme_light()

p2 <- ggplot(data, aes(x = QUAL)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(dp_summary$min_threshold_1st_QUAL, dp_summary$max_threshold_99th_QUAL), 
             color = "red", linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "QUAL score - raw data",
    x     = "DP (log10 scale)",
    y     = "Count"
  ) + theme_light()

png(filename = paste0(VIZ_DIR, "/01_0_filters.png"), width = 1000, height = 300, res = 100)
grid.arrange(p1, p2, nrow = 1)
dev.off()

rm(p1, p2) 

## FILTERING
data_filt <- data %>%
  filter(
    QUAL < dp_summary$max_threshold_99th_QUAL,
    QUAL > dp_summary$min_threshold_1st_QUAL,
    DP < dp_summary$max_threshold_99th,
    DP > dp_summary$min_threshold_1st
  )

# Save summary as as markdown table
dp_table <- rbind(
  "Raw"      = summary(data$DP),
  "Filtered" = summary(data_filt$DP)
)

qual_table <- rbind(
  "Raw"      = summary(data$QUAL),
  "Filtered" = summary(data_filt$QUAL)
)

kable(dp_table, format = "markdown", caption = "Read Depth (DP) Summary")
kable(qual_table, format = "markdown", caption = "PHRED Quality (QUAL) Summary")

rm(dp_table, qual_table, dp_summary)

# FIX order of chromosomes
unique(data_filt$CHROM)
data_filt <- data_filt %>%
  mutate(CHROM = factor(CHROM, levels = str_sort(unique(CHROM), numeric = TRUE)))

# ======================================================================================
p1 <- ggplot(data_filt, aes(x = QUAL)) +
  geom_histogram(aes(y = after_stat(count / max(count))), # normalizacia aby histogram bol do 1 ako CDF
                 bins = 50, fill = "gray85", colour = "white", linewidth = 0.2) +
  stat_ecdf(geom = "step", color = "palevioletred2", linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Normalized Count (Bars)",
    sec.axis = sec_axis(~., name = "Cumulative Probability (Line)")
  ) +
  labs(
    title = "QUAL Quality: Distribution & CDF",
    x = "QUAL (log10 scale)"
  ) + 
  theme_light()

p2 <- ggplot(data_filt, aes(x = DP)) +
  geom_histogram(aes(y = after_stat(count / max(count))), 
                 bins = 50, fill = "gray85", colour = "white", linewidth = 0.2) +
  stat_ecdf(geom = "step", color = "palevioletred2", linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Normalized Count (Bars)",
    sec.axis = sec_axis(~., name = "Cumulative Probability (Line)")
  ) +
  labs(
    title = "Read Depth (DP): Distribution & CDF",
    x = "DP (log10 scale)"
  ) + 
  theme_light()

p3 <- ggplot(data_filt, aes(x = QUAL, fill = TYPE, colour = TYPE)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  theme_light() +
  scale_fill_manual(values = type_colours) +
  scale_colour_manual(values = type_colours) +
  labs(
    title = "QUAL density – SNP vs INDEL",
    x     = "QUAL (log10 scale)",
    y     = "Density"
  )

p4 <- ggplot(data_filt, aes(x = TYPE, y = QUAL, fill = TYPE)) +
  geom_violin(outlier.alpha = 0.4, width = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4, fill = "white") +
  scale_y_log10() +
  scale_fill_manual(values = type_colours) +
  labs(
    title = "QUAL violin plot – SNP vs INDEL",
    x     = NULL,
    y     = "QUAL (log10 scale)"
  ) +
  theme_light()


p5 <- ggplot(data_filt, aes(x = DP, fill = TYPE, colour = TYPE)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  theme_light() +
  scale_fill_manual(values = type_colours) +
  scale_colour_manual(values = type_colours) +
  labs(
    title = "DP density – SNP vs INDEL",
    x     = "DP (log10 scale)",
    y     = "Density"
  )

p6 <- ggplot(data_filt, aes(x = TYPE, y = DP, fill = TYPE)) +
  geom_violin(outlier.alpha = 0.4, width = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4, fill = "white") +
  scale_y_log10() +
  scale_fill_manual(values = type_colours) +
  labs(
    title = "DP violin plot – SNP vs INDEL",
    x     = NULL,
    y     = "DP (log10 scale)"
  ) +
  theme_light()

my_layout <- rbind(
  c(1, 3, 4),
  c(2, 5, 6)
)
png(filename = paste0(VIZ_DIR, "/01_1_genome.png"), width = 1200, height = 600, res = 100)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = my_layout)
dev.off()

rm(p1, p2, p3, p4, p5, p6)

# ==============================================================================
# Exploratory analysis - by chromosome and variant type
summary_tbl <- data_filt %>%
  group_by(CHROM, TYPE) %>%
  summarise(
    n          = n(),
    QUAL_median = median(QUAL),
    QUAL_mean   = mean(QUAL),
    QUAL_IQR    = IQR(QUAL),
    DP_median   = median(DP),
    DP_mean     = mean(DP),
    DP_IQR      = IQR(DP),
    .groups     = "drop"
  )

p1 <- ggplot(summary_tbl, aes(x = TYPE, y = QUAL_mean, fill = TYPE)) +
  geom_violin(alpha = 0.5, outlier.shape = NA) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4, fill = "white") +
  geom_jitter(aes(color = TYPE), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10() +
  scale_fill_manual(values = type_colours) +
  scale_colour_manual(values = type_colours) +
  theme_light() +
  labs(
    title = "Spread of Median QUAL Scores Across Chromosomes",
    x     = "Variant Type",
    y     = "Median QUAL Score (log10 scale)"
  ) +
  theme(legend.position = "none") # Hides the legend since the x-axis already labels the types

p2 <- ggplot(summary_tbl, aes(x = CHROM, y = QUAL_median, color = TYPE)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = TYPE), alpha = 0.3) + 
  scale_y_log10() +
  scale_color_manual(values = type_colours) +
  scale_fill_manual(values = type_colours) +
  theme_light() +
  labs(
    title = "Median QUAL Score per Chromosome: SNP vs INDEL",
    x     = "Chromosome",
    y     = "Median QUAL Score (log10 scale)",
    color = "Variant Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


p3 <- ggplot(summary_tbl, aes(x = TYPE, y = DP_median, fill = TYPE)) +
  geom_violin(alpha = 0.5, outlier.shape = NA) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4, fill = "white") +
  geom_jitter(aes(color = TYPE), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10() +
  scale_fill_manual(values = type_colours) +
  scale_colour_manual(values = type_colours) +
  theme_light() +
  labs(
    title = "Spread of Median DP Scores Across Chromosomes",
    x     = "Variant Type",
    y     = "Median DP Score (log10 scale)"
  ) +
  theme(legend.position = "none") 


p4 <- ggplot(summary_tbl, aes(x = CHROM, y = DP_median, color = TYPE)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = TYPE), alpha = 0.3) + 
  scale_y_log10() +
  scale_color_manual(values = type_colours) +
  scale_fill_manual(values = type_colours) +
  theme_light() +
  labs(
    title = "Median DP Score per Chromosome: SNP vs INDEL",
    x     = "Chromosome",
    y     = "Median DP Score (log10 scale)",
    color = "Variant Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


png(filename = paste0(VIZ_DIR, "/01_2_chr.png"), width = 1000, height = 600, res = 100)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

rm(p1, p2, p3, p4, summary_tbl)






# =============================================================================
# Task 9

dim(data_filt)

cor_overall_s <- cor.test(log10(data_filt$QUAL), log10(data_filt$DP),
                        method = "spearman")

cor_overall_p <- cor.test(log10(data_filt$QUAL), log10(data_filt$DP),
                          method = "pearson")

cor_by_type <- data_filt %>%
  group_by(TYPE) %>%
  summarise(
    spearman_rho = cor(log10(QUAL), log10(DP), method = "spearman"),
    spearman_p   = cor.test(log10(QUAL), log10(DP), method = "spearman", exact = FALSE)$p.value,
    pearson_r    = cor(log10(QUAL), log10(DP), method = "pearson"),
    pearson_p    = cor.test(log10(QUAL), log10(DP), method = "pearson")$p.value,
    
    .groups = "drop",
    n = n()
  )



print(cor_overall_s)
print(cor_overall_p)
print(cor_by_type)



p <- ggplot(data_filt, aes(x = DP, y = QUAL, color = TYPE)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_y_log10() + 
  scale_x_log10() +
  
  scale_color_manual(values = type_colours) +
  scale_fill_manual(values = type_colours) +
  
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
  
  theme_light() +
  theme(legend.position = "bottom") + 
  labs(
    title = "Read Depth vs. PHRED Quality",
    x = "Read Depth (DP) (log10 scale)",
    y = "QUAL (log10 scale)",
    color = "Variant Type"
  )

p_marginal <- ggMarginal(
  p, 
  type = "density", 
  groupColour = TRUE, 
  groupFill = TRUE, 
  alpha = 0.4
)

png(filename = paste0(VIZ_DIR, "/01_3_scatter.png"), width = 600, height = 400, res = 100)
print(p_marginal)
dev.off()



p_sub <- ggplot(data_filt, aes(x = DP, y = QUAL, colour = TYPE)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_manual(values = type_colours) +
  scale_fill_manual(values = type_colours) +
  
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
  
  facet_wrap(~ TYPE) +
  
  theme_light() +
  labs(
    title = "Read Depth vs. QUAL score",
    x = "Read Depth (log10 scale)",
    y = "QUAL Score (log10 scale)"
  )

png(filename = paste0(VIZ_DIR, "/01_4_scatter_sub.png"), width = 800, height = 400, res = 100)
print(p_sub)
dev.off()
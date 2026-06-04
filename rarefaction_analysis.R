library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)

# ---------------------------------------------------------
# RAREFACTION CURVES: BASAL VS FINAL (ALL SAMPLES)
# WITH OBSERVED AND CHAO1 RICHNESS ESTIMATES
# ---------------------------------------------------------

# ---------------------------------------------------------
# 1. Subset phyloseq object by timepoint
# ---------------------------------------------------------
ps_basal <- subset_samples(ps, timepoint == "Basal")
ps_final <- subset_samples(ps, timepoint == "Final")

# ---------------------------------------------------------
# 2. Function to compute mean rarefaction curve per timepoint
# ---------------------------------------------------------
compute_mean_rarefaction <- function(ps_obj, timepoint_label) {
  
  otu <- as.data.frame(t(otu_table(ps_obj)))
  all_individual <- NULL
  
  for (sample_name in rownames(otu)) {
    sample_reads <- as.numeric(otu[sample_name, ])
    sample_total <- sum(sample_reads)
    
    if (sample_total > 0) {
      read_depths <- seq(100, sample_total, by = 100)
      richness_vals <- sapply(read_depths, function(x) {
        if (x <= sample_total) {
          rarefy(sample_reads, x)
        } else {
          NA
        }
      })
      
      sample_df <- tibble(
        SampleID = sample_name,
        Reads = read_depths,
        Richness = as.numeric(richness_vals),
        timepoint = timepoint_label,
        line_type = "individual"
      )
      
      all_individual <- bind_rows(all_individual, sample_df)
    }
  }
  
  # Compute MEAN rarefaction curve for this timepoint
  mean_rare <- all_individual %>%
    group_by(Reads) %>%
    summarise(
      Richness = mean(Richness, na.rm = TRUE),
      timepoint = first(timepoint),
      line_type = "mean",
      .groups = "drop"
    )
  
  combined <- bind_rows(all_individual, mean_rare)
  return(combined)
}

# ---------------------------------------------------------
# 3. Generate rarefaction data for both timepoints
# ---------------------------------------------------------
cat("Computing rarefaction curves...\n")
rare_basal <- compute_mean_rarefaction(ps_basal, "Basal")
rare_final <- compute_mean_rarefaction(ps_final, "Final")
rare_all <- bind_rows(rare_basal, rare_final)

cat("Rarefaction data computed successfully.\n")

# ---------------------------------------------------------
# 4. Compute Observed and Chao1 richness for all samples
# ---------------------------------------------------------
cat("Computing Observed and Chao1 richness estimates...\n")

richness_df <- estimate_richness(ps, measures = c("Observed", "Chao1")) %>%
  rownames_to_column("SampleID")

meta <- data.frame(sample_data(ps)) %>%
  rownames_to_column("SampleID")

richness_df <- richness_df %>%
  left_join(meta, by = "SampleID")

# ---------------------------------------------------------
# 5. Timepoint colors - dark primary, light secondary
# ---------------------------------------------------------
timepoint_colors <- c(
  Basal = "#C0392B",    # dark red
  Final = "#2874A6"     # dark blue
)

# ---------------------------------------------------------
# 6. PLOT 1: Rarefaction curves - BASAL vs FINAL
# ---------------------------------------------------------
p_rare_comparison <- ggplot() +
  # Individual sample curves (light, desaturated)
  geom_line(data = filter(rare_all, line_type == "individual"), 
            aes(x = Reads, y = Richness, color = timepoint, group = SampleID),
            linewidth = 0.5, alpha = 0.4) +
  # Mean curves (dark, bold)
  geom_line(data = filter(rare_all, line_type == "mean"),
            aes(x = Reads, y = Richness, color = timepoint, linetype = timepoint),
            linewidth = 1.5, alpha = 1) +
  scale_color_manual(values = timepoint_colors) +
  scale_linetype_manual(values = c("Basal" = "solid", "Final" = "dashed")) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "Rarefaction Curves: Basal vs Final Timepoint\n(Individual samples light, Mean curves bold)",
    x = "Number of Reads Sampled",
    y = "Observed OTU Richness",
    color = "Timepoint",
    linetype = "Timepoint"
  )

print(p_rare_comparison)

# ---------------------------------------------------------
# 7. PLOT 2: Observed Richness by Timepoint
# ---------------------------------------------------------
p_observed <- ggplot(richness_df, aes(x = timepoint, y = Observed, 
                                       fill = timepoint)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "Observed Richness: Basal vs Final",
    x = "Timepoint",
    y = "Observed OTU Count",
    fill = "Timepoint"
  )

print(p_observed)

# ---------------------------------------------------------
# 8. PLOT 3: Chao1 Richness by Timepoint
# ---------------------------------------------------------
p_chao1 <- ggplot(richness_df, aes(x = timepoint, y = Chao1, 
                                    fill = timepoint)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "Chao1 Richness: Basal vs Final",
    x = "Timepoint",
    y = "Chao1 Richness Estimate",
    fill = "Timepoint"
  )

print(p_chao1)

# ---------------------------------------------------------
# 9. SAVE PLOTS
# ---------------------------------------------------------
outdir <- "C:/Users/mgorton/OneDrive - Oklahoma A and M System/Desktop/16s rerun/new results"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

cat("Saving plots to:", outdir, "\n")

ggsave(file.path(outdir, "Rarefaction_Basal_vs_Final.png"),
       plot = p_rare_comparison, width = 11, height = 7, dpi = 300)

ggsave(file.path(outdir, "Observed_Richness_Basal_vs_Final.png"),
       plot = p_observed, width = 8, height = 6, dpi = 300)

ggsave(file.path(outdir, "Chao1_Richness_Basal_vs_Final.png"),
       plot = p_chao1, width = 8, height = 6, dpi = 300)

cat("Plots saved successfully.\n")

# ---------------------------------------------------------
# 10. SUMMARY STATISTICS
# ---------------------------------------------------------
rare_summary <- rare_all %>%
  filter(line_type == "mean") %>%
  group_by(timepoint) %>%
  summarise(
    max_reads = max(Reads, na.rm = TRUE),
    final_richness = last(Richness),
    .groups = "drop"
  )

cat("\n=== Rarefaction Summary (Mean Curves) ===\n")
print(rare_summary)

# Observed richness summary
obs_summary <- richness_df %>%
  group_by(timepoint) %>%
  summarise(
    n = n(),
    mean_observed = mean(Observed, na.rm = TRUE),
    sd_observed = sd(Observed, na.rm = TRUE),
    median_observed = median(Observed, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Observed Richness Summary ===\n")
print(obs_summary)

# Chao1 richness summary
chao_summary <- richness_df %>%
  group_by(timepoint) %>%
  summarise(
    n = n(),
    mean_chao1 = mean(Chao1, na.rm = TRUE),
    sd_chao1 = sd(Chao1, na.rm = TRUE),
    median_chao1 = median(Chao1, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Chao1 Richness Summary ===\n")
print(chao_summary)

# Save all summaries
write.csv(rare_summary,
          file.path(outdir, "Rarefaction_Summary_Basal_vs_Final.csv"),
          row.names = FALSE)

write.csv(obs_summary,
          file.path(outdir, "Observed_Richness_Summary.csv"),
          row.names = FALSE)

write.csv(chao_summary,
          file.path(outdir, "Chao1_Richness_Summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 11. STATISTICAL TESTS: BASAL vs FINAL
# ---------------------------------------------------------

cat("\n=========================================\n")
cat("STATISTICAL COMPARISON: BASAL vs FINAL\n")
cat("=========================================\n")

basal_obs <- richness_df %>% filter(timepoint == "Basal") %>% pull(Observed)
final_obs <- richness_df %>% filter(timepoint == "Final") %>% pull(Observed)

basal_chao <- richness_df %>% filter(timepoint == "Basal") %>% pull(Chao1)
final_chao <- richness_df %>% filter(timepoint == "Final") %>% pull(Chao1)

# Mann-Whitney U tests
mw_obs <- wilcox.test(basal_obs, final_obs, paired = FALSE)
mw_chao <- wilcox.test(basal_chao, final_chao, paired = FALSE)

cat("\n--- Observed Richness (Mann-Whitney U) ---\n")
cat("Basal mean:", round(mean(basal_obs, na.rm = TRUE), 2), "\n")
cat("Final mean:", round(mean(final_obs, na.rm = TRUE), 2), "\n")
cat("P-value:", round(mw_obs$p.value, 4), "\n")
if (mw_obs$p.value < 0.05) {
  cat("Result: SIGNIFICANT difference (p < 0.05)\n")
} else {
  cat("Result: NO significant difference (p >= 0.05)\n")
}

cat("\n--- Chao1 Richness (Mann-Whitney U) ---\n")
cat("Basal mean:", round(mean(basal_chao, na.rm = TRUE), 2), "\n")
cat("Final mean:", round(mean(final_chao, na.rm = TRUE), 2), "\n")
cat("P-value:", round(mw_chao$p.value, 4), "\n")
if (mw_chao$p.value < 0.05) {
  cat("Result: SIGNIFICANT difference (p < 0.05)\n")
} else {
  cat("Result: NO significant difference (p >= 0.05)\n")
}

cat("\n=========================================\n")
cat("✓ Rarefaction & Richness analysis complete!\n")
cat("✓ All plots and summaries saved to:\n")
cat("  ", outdir, "\n")
cat("=========================================\n")


# ===================================================================================
# SECTION 2: ALPHA DIVERSITY (SHANNON, SIMPSON, ACE) - BASAL VS FINAL
# ===================================================================================

cat("\n\n")
cat("===================================================================================\n")
cat("SECTION 2: COMPREHENSIVE ALPHA DIVERSITY ANALYSIS\n")
cat("===================================================================================\n\n")

library(FSA)
library(dplyr)
library(openxlsx)

# ---------------------------------------------------------
# 1. Compute all alpha diversity metrics
# ---------------------------------------------------------
cat("Computing Shannon, Simpson, and ACE diversity estimates...\n")

alpha_metrics <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "ACE", "Chao1")) %>%
  rownames_to_column("SampleID")

alpha_df <- alpha_metrics %>%
  left_join(meta, by = "SampleID")

# ---------------------------------------------------------
# 2. PLOT 4: Shannon Diversity - BASAL vs FINAL
# ---------------------------------------------------------
p_shannon <- ggplot(alpha_df, aes(x = timepoint, y = Shannon, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "Shannon Diversity: Basal vs Final",
    x = "Timepoint",
    y = "Shannon Diversity Index",
    fill = "Timepoint"
  )

print(p_shannon)

# ---------------------------------------------------------
# 3. PLOT 5: Simpson Diversity - BASAL vs FINAL
# ---------------------------------------------------------
p_simpson <- ggplot(alpha_df, aes(x = timepoint, y = Simpson, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "Simpson Diversity: Basal vs Final",
    x = "Timepoint",
    y = "Simpson Diversity Index",
    fill = "Timepoint"
  )

print(p_simpson)

# ---------------------------------------------------------
# 4. PLOT 6: ACE Richness - BASAL vs FINAL
# ---------------------------------------------------------
p_ace <- ggplot(alpha_df, aes(x = timepoint, y = ACE, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "ACE Richness: Basal vs Final",
    x = "Timepoint",
    y = "ACE Richness Estimate",
    fill = "Timepoint"
  )

print(p_ace)

# ---------------------------------------------------------
# 5. Save alpha diversity plots
# ---------------------------------------------------------
ggsave(file.path(outdir, "Shannon_Diversity_Basal_vs_Final.png"),
       plot = p_shannon, width = 8, height = 6, dpi = 300)

ggsave(file.path(outdir, "Simpson_Diversity_Basal_vs_Final.png"),
       plot = p_simpson, width = 8, height = 6, dpi = 300)

ggsave(file.path(outdir, "ACE_Richness_Basal_vs_Final.png"),
       plot = p_ace, width = 8, height = 6, dpi = 300)

cat("✓ Alpha diversity plots saved.\n")

# ---------------------------------------------------------
# 6. Statistical tests for all alpha metrics
# ---------------------------------------------------------
cat("\nRunning statistical tests for alpha diversity metrics...\n")

basal_shannon <- alpha_df %>% filter(timepoint == "Basal") %>% pull(Shannon)
final_shannon <- alpha_df %>% filter(timepoint == "Final") %>% pull(Shannon)

basal_simpson <- alpha_df %>% filter(timepoint == "Basal") %>% pull(Simpson)
final_simpson <- alpha_df %>% filter(timepoint == "Final") %>% pull(Simpson)

basal_ace <- alpha_df %>% filter(timepoint == "Basal") %>% pull(ACE)
final_ace <- alpha_df %>% filter(timepoint == "Final") %>% pull(ACE)

# Mann-Whitney U tests for all metrics
mw_shannon <- wilcox.test(basal_shannon, final_shannon, paired = FALSE)
mw_simpson <- wilcox.test(basal_simpson, final_simpson, paired = FALSE)
mw_ace <- wilcox.test(basal_ace, final_ace, paired = FALSE)

# Create summary table
alpha_stats <- data.frame(
  Metric = c("Observed", "Chao1", "Shannon", "Simpson", "ACE"),
  Basal_Mean = c(
    mean(basal_obs, na.rm = TRUE),
    mean(basal_chao, na.rm = TRUE),
    mean(basal_shannon, na.rm = TRUE),
    mean(basal_simpson, na.rm = TRUE),
    mean(basal_ace, na.rm = TRUE)
  ),
  Final_Mean = c(
    mean(final_obs, na.rm = TRUE),
    mean(final_chao, na.rm = TRUE),
    mean(final_shannon, na.rm = TRUE),
    mean(final_simpson, na.rm = TRUE),
    mean(final_ace, na.rm = TRUE)
  ),
  Basal_SD = c(
    sd(basal_obs, na.rm = TRUE),
    sd(basal_chao, na.rm = TRUE),
    sd(basal_shannon, na.rm = TRUE),
    sd(basal_simpson, na.rm = TRUE),
    sd(basal_ace, na.rm = TRUE)
  ),
  Final_SD = c(
    sd(final_obs, na.rm = TRUE),
    sd(final_chao, na.rm = TRUE),
    sd(final_shannon, na.rm = TRUE),
    sd(final_simpson, na.rm = TRUE),
    sd(final_ace, na.rm = TRUE)
  ),
  P_Value = c(
    mw_obs$p.value,
    mw_chao$p.value,
    mw_shannon$p.value,
    mw_simpson$p.value,
    mw_ace$p.value
  ),
  Significant = c(
    ifelse(mw_obs$p.value < 0.05, "Yes", "No"),
    ifelse(mw_chao$p.value < 0.05, "Yes", "No"),
    ifelse(mw_shannon$p.value < 0.05, "Yes", "No"),
    ifelse(mw_simpson$p.value < 0.05, "Yes", "No"),
    ifelse(mw_ace$p.value < 0.05, "Yes", "No")
  )
)

cat("\n=== ALPHA DIVERSITY STATISTICS ===\n")
print(alpha_stats)

# ---------------------------------------------------------
# 7. BETA DIVERSITY: PCoA Analysis
# ---------------------------------------------------------
cat("\n\nComputing beta diversity metrics...\n")

dist_bray <- phyloseq::distance(ps, method = "bray")
dist_unifrac <- phyloseq::distance(ps, method = "unifrac", weighted = FALSE)
dist_wunifrac <- phyloseq::distance(ps, method = "unifrac", weighted = TRUE)

ord_bray <- ordinate(ps, method = "PCoA", distance = dist_bray)
ord_unifrac <- ordinate(ps, method = "PCoA", distance = dist_unifrac)
ord_wunifrac <- ordinate(ps, method = "PCoA", distance = dist_wunifrac)

# PERMANOVA tests
meta_df <- data.frame(sample_data(ps))

perm_bray <- adonis2(dist_bray ~ timepoint, data = meta_df)
perm_unifrac <- adonis2(dist_unifrac ~ timepoint, data = meta_df)
perm_wunifrac <- adonis2(dist_wunifrac ~ timepoint, data = meta_df)

p_bray_perm <- round(perm_bray$`Pr(>F)`[1], 4)
p_unifrac_perm <- round(perm_unifrac$`Pr(>F)`[1], 4)
p_wunifrac_perm <- round(perm_wunifrac$`Pr(>F)`[1], 4)

# ---------------------------------------------------------
# 8. PLOT 7: PCoA - Bray-Curtis
# ---------------------------------------------------------
p_bray_ord <- plot_ordination(ps, ord_bray, color = "timepoint") +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(aes(color = timepoint), type = "t", linewidth = 1.2) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(title = paste0("PCoA – Bray-Curtis\n(PERMANOVA p = ", p_bray_perm, ")"),
       color = "Timepoint")

print(p_bray_ord)

# ---------------------------------------------------------
# 9. PLOT 8: PCoA - Unweighted UniFrac
# ---------------------------------------------------------
p_unifrac_ord <- plot_ordination(ps, ord_unifrac, color = "timepoint") +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(aes(color = timepoint), type = "t", linewidth = 1.2) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(title = paste0("PCoA – Unweighted UniFrac\n(PERMANOVA p = ", p_unifrac_perm, ")"),
       color = "Timepoint")

print(p_unifrac_ord)

# ---------------------------------------------------------
# 10. PLOT 9: PCoA - Weighted UniFrac
# ---------------------------------------------------------
p_wunifrac_ord <- plot_ordination(ps, ord_wunifrac, color = "timepoint") +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(aes(color = timepoint), type = "t", linewidth = 1.2) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(title = paste0("PCoA – Weighted UniFrac\n(PERMANOVA p = ", p_wunifrac_perm, ")"),
       color = "Timepoint")

print(p_wunifrac_ord)

# ---------------------------------------------------------
# 11. Save beta diversity plots
# ---------------------------------------------------------
ggsave(file.path(outdir, "PCoA_BrayCurtis_Basal_vs_Final.png"),
       plot = p_bray_ord, width = 9, height = 7, dpi = 300)

ggsave(file.path(outdir, "PCoA_UnweightedUniFrac_Basal_vs_Final.png"),
       plot = p_unifrac_ord, width = 9, height = 7, dpi = 300)

ggsave(file.path(outdir, "PCoA_WeightedUniFrac_Basal_vs_Final.png"),
       plot = p_wunifrac_ord, width = 9, height = 7, dpi = 300)

cat("✓ Beta diversity plots saved.\n")

# ---------------------------------------------------------
# 12. Create beta diversity statistics table
# ---------------------------------------------------------
beta_stats <- data.frame(
  Distance_Metric = c("Bray-Curtis", "Unweighted UniFrac", "Weighted UniFrac"),
  PERMANOVA_F = c(perm_bray$F[1], perm_unifrac$F[1], perm_wunifrac$F[1]),
  PERMANOVA_P = c(perm_bray$`Pr(>F)`[1], perm_unifrac$`Pr(>F)`[1], perm_wunifrac$`Pr(>F)`[1]),
  Significant = c(
    ifelse(perm_bray$`Pr(>F)`[1] < 0.05, "Yes", "No"),
    ifelse(perm_unifrac$`Pr(>F)`[1] < 0.05, "Yes", "No"),
    ifelse(perm_wunifrac$`Pr(>F)`[1] < 0.05, "Yes", "No")
  )
)

cat("\n=== BETA DIVERSITY STATISTICS (PERMANOVA) ===\n")
print(beta_stats)

# ---------------------------------------------------------
# 13. Save all statistics to Excel
# ---------------------------------------------------------
cat("\nSaving all statistics to Excel...\n")

outfile_stats <- file.path(outdir, "Alpha_Beta_Diversity_Statistics_Basal_vs_Final.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "Alpha_Diversity")
writeData(wb, "Alpha_Diversity", alpha_stats)

addWorksheet(wb, "Beta_Diversity")
writeData(wb, "Beta_Diversity", beta_stats)

addWorksheet(wb, "Raw_Data")
writeData(wb, "Raw_Data", alpha_df)

saveWorkbook(wb, outfile_stats, overwrite = TRUE)

cat("✓ Statistics exported to Excel.\n")

cat("\n\n")
cat("===================================================================================\n")
cat("✓ SECTION 2: COMPREHENSIVE ALPHA/BETA DIVERSITY ANALYSIS COMPLETE!\n")
cat("✓ All plots and statistics saved to:\n")
cat("  ", outdir, "\n")
cat("===================================================================================\n")

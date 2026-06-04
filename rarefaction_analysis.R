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
  
  # Get OTU table (samples x taxa)
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
        line_type = "individual"  # Mark as individual sample
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
  
  # Combine individual + mean
  combined <- bind_rows(all_individual, mean_rare)
  
  return(combined)
}

# ---------------------------------------------------------
# 3. Generate rarefaction data for both timepoints
# ---------------------------------------------------------
cat("Computing rarefaction curves for Basal samples...\n")
rare_basal <- compute_mean_rarefaction(ps_basal, "Basal")

cat("Computing rarefaction curves for Final samples...\n")
rare_final <- compute_mean_rarefaction(ps_final, "Final")

# Combine both
rare_all <- bind_rows(rare_basal, rare_final)

cat("Rarefaction data computed successfully.\n")

# ---------------------------------------------------------
# 4. Compute Observed and Chao1 richness for all samples
# ---------------------------------------------------------
cat("Computing Observed and Chao1 richness estimates...\n")

richness_df <- estimate_richness(ps, measures = c("Observed", "Chao1")) %>%
  rownames_to_column("SampleID")

# Add metadata
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

timepoint_colors_light <- c(
  Basal = "#F5B7B1",    # light red (saturated)
  Final = "#AED6F1"     # light blue (saturated)
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

cat("All plots saved successfully.\n")

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
cat("✓ All analyses complete!\n")
cat("✓ Plots and summaries saved to:\n")
cat("  ", outdir, "\n")
cat("=========================================\n")

## =============================================================================
##  Transcriptome-AutoPipeline — Statistical analysis of snoRNA-containing
##  retained-intron (RI) events  (replication-aware version)
##
##  Input : results/RI_events_with_snoRNA.tsv   (SUPPA2 PSI, events x SRR runs)
##          metadata/SraRunTable.csv            (SRA run table, PRJNA751490)
##  Output: results/stats/*.tsv  and  figures/*.png
##
##  Design decisions (fixed BEFORE looking at results):
##   * Only bulk libraries are used (LibrarySource == "TRANSCRIPTOMIC");
##     the 12 single-cell runs are a different assay and are excluded.
##   * Every library (BioSample) was sequenced on 4 lanes = 4 SRR runs.
##     Lanes are technical replicates, so PSI is averaged per BioSample.
##     The statistical unit is the library, not the run.
##   * Event rows duplicated only because one intron hosts two snoRNAs are
##     tested once (snoRNA names are joined, e.g. "SNORD49A;SNORD49B").
##   * Groups: Healthy / Normal (adjacent) / Low / High, where Low/High is the
##     authors' `grade` annotation (low_grade = Gleason 6-7, high_grade = >= 8).
##   * Patient ID is not in the run table; a proxy (age + Gleason string)
##     links tumour and adjacent-normal libraries of the same patient.
##   * Multiple testing: Benjamini-Hochberg within each analysis and contrast.
##   * All seven analyses below are reported, significant or not.
##
##  Analyses
##   S1  Per-event Wilcoxon rank-sum, all events, 4 contrasts
##   S2  Same as S1 after label-free independent filtering
##   S3  snoRNA-level Wilcoxon (PSI averaged over the events of each snoRNA)
##   S4  Spearman correlation of PSI with Gleason score (tumour libraries)
##   S5  Per-event linear mixed model on lane-level data, PSI ~ group + (1|BioSample)
##   S6  Global RI index (mean z-scored PSI per library): 4 contrasts + paired T vs N
##   S7  Paired Wilcoxon signed-rank, tumour vs adjacent normal, per event
##
##  Required packages: readr dplyr tidyr stringr purrr tibble ggplot2 lme4 lmerTest
## =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(ggplot2); library(lmerTest)
})

args      <- commandArgs(trailingOnly = TRUE)
PSI_FILE  <- if (length(args) >= 1) args[1] else "results/RI_events_with_snoRNA.tsv"
META_FILE <- if (length(args) >= 2) args[2] else "metadata/SraRunTable.csv"
OUT_DIR   <- "results/stats"
FIG_DIR   <- "figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_PRESENT <- 0.80   # S2 filter: event measured in >= 80% of libraries
MIN_SD      <- 0.01   # S2 filter: SD of PSI across libraries >= 0.01
CONTRASTS <- tibble(
  contrast = c("Low_vs_Healthy", "High_vs_Healthy", "Normal_vs_Healthy", "Low_vs_High"),
  g1       = c("Low",  "High",    "Normal",  "Low"),
  g2       = c("Healthy", "Healthy", "Healthy", "High")
)

## ── 1. Metadata: bulk only, groups, patient proxy ─────────────────────────────
meta <- read_csv(META_FILE, show_col_types = FALSE) %>%
  filter(LibrarySource == "TRANSCRIPTOMIC") %>%
  transmute(
    Run, BioSample, condition, gleason, age,
    group = case_when(
      str_detect(tolower(condition), "healthy") ~ "Healthy",
      str_detect(tolower(condition), "normal")  ~ "Normal",
      str_detect(tolower(grade), "low")         ~ "Low",
      str_detect(tolower(grade), "high")        ~ "High",
      TRUE ~ NA_character_
    ),
    gleason_score = suppressWarnings(as.integer(
      ifelse(str_detect(gleason, "="),
             str_extract(gleason, "(?<==)\\s*\\d+"),
             map_chr(str_extract_all(gleason, "\\d+"), ~ as.character(sum(as.integer(.x)))))
    )),
    patient = ifelse(is.na(gleason), BioSample, paste(age, gleason, sep = "|"))
  ) %>%
  filter(!is.na(group))

## ── 2. PSI: unique events, long format, lane-level and library-level ──────────
psi_raw <- read_tsv(PSI_FILE, show_col_types = FALSE)
snorna_map <- psi_raw %>% distinct(eventID, snoRNAName)
event_names <- snorna_map %>% group_by(eventID) %>%
  summarise(snoRNAName = paste(sort(unique(snoRNAName)), collapse = ";"), .groups = "drop")

run_long <- psi_raw %>%
  distinct(eventID, .keep_all = TRUE) %>%
  select(-snoRNAName) %>%
  pivot_longer(-eventID, names_to = "Run", values_to = "PSI") %>%
  inner_join(meta, by = "Run")

lib_long <- run_long %>%
  group_by(eventID, BioSample, group, condition, patient, gleason_score) %>%
  summarise(PSI = if (all(is.na(PSI))) NA_real_ else mean(PSI, na.rm = TRUE), .groups = "drop")

libs <- lib_long %>% distinct(BioSample, group, condition, patient)
cat("Bulk runs:", n_distinct(run_long$Run), "| libraries:", nrow(libs),
    "| unique events:", n_distinct(lib_long$eventID), "\n")
print(table(libs$group))

## ── helpers ──────────────────────────────────────────────────────────────────
wilcox_two <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3)
    return(tibble(n1 = length(x), n2 = length(y), deltaPSI = NA_real_,
                  HL_estimate = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_))
  w <- tryCatch(suppressWarnings(wilcox.test(x, y, conf.int = TRUE, exact = FALSE)),
                error = function(e) suppressWarnings(wilcox.test(x, y, exact = FALSE)))
  tibble(n1 = length(x), n2 = length(y), deltaPSI = mean(x) - mean(y),
         HL_estimate = if (!is.null(w$estimate)) unname(w$estimate) else NA_real_,
         CI_low  = if (!is.null(w$conf.int)) w$conf.int[1] else NA_real_,
         CI_high = if (!is.null(w$conf.int)) w$conf.int[2] else NA_real_,
         p = w$p.value)
}

per_feature_contrasts <- function(df, feature) {
  pmap_dfr(CONTRASTS, function(contrast, g1, g2) {
    df %>% group_by(.data[[feature]]) %>%
      group_modify(~ wilcox_two(.x$PSI[.x$group == g1], .x$PSI[.x$group == g2])) %>%
      ungroup() %>% mutate(contrast = contrast, FDR = p.adjust(p, "BH"), .before = 1)
  })
}

sig_summary <- function(res, analysis) {
  res %>% group_by(contrast) %>%
    summarise(analysis = analysis, n_tested = sum(!is.na(p)),
              n_p_lt_0.05 = sum(p < 0.05, na.rm = TRUE),
              n_FDR_lt_0.05 = sum(FDR < 0.05, na.rm = TRUE),
              min_FDR = suppressWarnings(min(FDR, na.rm = TRUE)), .groups = "drop") %>%
    relocate(analysis)
}

## ── S1. Per-event Wilcoxon, all events ───────────────────────────────────────
s1 <- per_feature_contrasts(lib_long, "eventID") %>% left_join(event_names, by = "eventID")
write_tsv(s1, file.path(OUT_DIR, "S1_event_wilcoxon_all.tsv"))

## ── S2. Independent (label-free) filtering, then per-event Wilcoxon ──────────
keep <- lib_long %>% group_by(eventID) %>%
  summarise(present = mean(!is.na(PSI)), sd = sd(PSI, na.rm = TRUE), .groups = "drop") %>%
  filter(present >= MIN_PRESENT, sd >= MIN_SD) %>% pull(eventID)
cat("S2 filter kept", length(keep), "of", n_distinct(lib_long$eventID), "events\n")
s2 <- per_feature_contrasts(filter(lib_long, eventID %in% keep), "eventID") %>%
  left_join(event_names, by = "eventID")
write_tsv(s2, file.path(OUT_DIR, "S2_event_wilcoxon_filtered.tsv"))

## ── S3. snoRNA-level aggregation ─────────────────────────────────────────────
sno_long <- lib_long %>% inner_join(snorna_map, by = "eventID", relationship = "many-to-many") %>%
  group_by(snoRNAName, BioSample, group) %>%
  summarise(PSI = if (all(is.na(PSI))) NA_real_ else mean(PSI, na.rm = TRUE), .groups = "drop")
s3 <- per_feature_contrasts(sno_long, "snoRNAName")
write_tsv(s3, file.path(OUT_DIR, "S3_snoRNA_wilcoxon.tsv"))

## ── S4. Spearman with Gleason score (tumour libraries only) ──────────────────
s4 <- lib_long %>% filter(condition == "Tumor", !is.na(gleason_score), !is.na(PSI)) %>%
  group_by(eventID) %>% filter(n() >= 10, sd(PSI) > 0) %>%
  summarise(n = n(),
            rho = suppressWarnings(cor(PSI, gleason_score, method = "spearman")),
            p = suppressWarnings(cor.test(PSI, gleason_score, method = "spearman", exact = FALSE)$p.value),
            .groups = "drop") %>%
  mutate(FDR = p.adjust(p, "BH"), contrast = "Spearman_vs_Gleason") %>%
  left_join(event_names, by = "eventID") %>% arrange(p)
write_tsv(s4, file.path(OUT_DIR, "S4_spearman_gleason_tumour.tsv"))

## ── S5. Linear mixed model on lane-level data ────────────────────────────────
fit_lmm <- function(d, g1, g2) {
  d <- d %>% filter(group %in% c(g1, g2), !is.na(PSI)) %>%
    mutate(group = factor(group, levels = c(g2, g1)))
  if (n_distinct(d$group) < 2 || n_distinct(d$BioSample) < 6 || sd(d$PSI) == 0)
    return(tibble(estimate = NA_real_, p = NA_real_))
  m <- tryCatch(suppressMessages(suppressWarnings(lmer(PSI ~ group + (1 | BioSample), data = d))),
                error = function(e) NULL)
  if (is.null(m)) return(tibble(estimate = NA_real_, p = NA_real_))
  co <- summary(m)$coefficients
  tibble(estimate = co[2, "Estimate"], p = co[2, "Pr(>|t|)"])
}
s5 <- pmap_dfr(CONTRASTS, function(contrast, g1, g2) {
  run_long %>% group_by(eventID) %>% group_modify(~ fit_lmm(.x, g1, g2)) %>% ungroup() %>%
    mutate(contrast = contrast, FDR = p.adjust(p, "BH"), .before = 1)
}) %>% left_join(event_names, by = "eventID")
write_tsv(s5, file.path(OUT_DIR, "S5_event_mixed_model.tsv"))

## ── S6. Global RI index ──────────────────────────────────────────────────────
ri_index <- lib_long %>% filter(eventID %in% keep) %>%
  group_by(eventID) %>% mutate(z = (PSI - mean(PSI, na.rm = TRUE)) / sd(PSI, na.rm = TRUE)) %>%
  ungroup() %>% group_by(BioSample, group, condition, patient) %>%
  summarise(RI_index = mean(z, na.rm = TRUE), .groups = "drop")
write_tsv(ri_index, file.path(OUT_DIR, "S6_RI_index_per_library.tsv"))

s6_unpaired <- pmap_dfr(CONTRASTS, function(contrast, g1, g2)
  wilcox_two(ri_index$RI_index[ri_index$group == g1], ri_index$RI_index[ri_index$group == g2]) %>%
    mutate(contrast = contrast, .before = 1))
pairs_idx <- ri_index %>% filter(condition %in% c("Tumor", "Normal")) %>%
  group_by(patient, condition) %>% summarise(RI_index = mean(RI_index), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = RI_index) %>% drop_na(Tumor, Normal)
w6 <- wilcox.test(pairs_idx$Tumor, pairs_idx$Normal, paired = TRUE, exact = FALSE)
s6 <- bind_rows(s6_unpaired,
                tibble(contrast = "Tumor_vs_Normal_paired", n1 = nrow(pairs_idx), n2 = nrow(pairs_idx),
                       deltaPSI = median(pairs_idx$Tumor - pairs_idx$Normal), p = w6$p.value)) %>%
  rename(delta_index = deltaPSI) %>% mutate(FDR = p.adjust(p, "BH"))
write_tsv(s6, file.path(OUT_DIR, "S6_RI_index_tests.tsv"))

## ── S7. Paired tumour vs adjacent normal, per event ──────────────────────────
pairs_ev <- lib_long %>% filter(condition %in% c("Tumor", "Normal")) %>%
  group_by(eventID, patient, condition) %>% summarise(PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = PSI) %>%
  filter(is.finite(Tumor), is.finite(Normal))
s7 <- pairs_ev %>% group_by(eventID) %>%
  filter(n() >= 6, any(Tumor != Normal)) %>%
  summarise(n_pairs = n(), median_dPSI_T_minus_N = median(Tumor - Normal),
            p = suppressWarnings(wilcox.test(Tumor, Normal, paired = TRUE, exact = FALSE)$p.value),
            .groups = "drop") %>%
  mutate(contrast = "Tumor_vs_Normal_paired", FDR = p.adjust(p, "BH")) %>%
  left_join(event_names, by = "eventID") %>% arrange(p)
write_tsv(s7, file.path(OUT_DIR, "S7_paired_tumour_vs_normal.tsv"))

## ── Summary table ────────────────────────────────────────────────────────────
summary_tbl <- bind_rows(
  sig_summary(s1, "S1 per-event Wilcoxon (all events)"),
  sig_summary(s2, "S2 per-event Wilcoxon (filtered)"),
  sig_summary(s3, "S3 snoRNA-level Wilcoxon"),
  sig_summary(s4, "S4 Spearman vs Gleason (tumour)"),
  sig_summary(s5, "S5 mixed model (lane-level)"),
  sig_summary(s6, "S6 global RI index"),
  sig_summary(s7, "S7 paired tumour vs normal")
)
write_tsv(summary_tbl, file.path(OUT_DIR, "SUMMARY_all_analyses.tsv"))
cat("\n===== SUMMARY =====\n"); print(as.data.frame(summary_tbl), row.names = FALSE)

cat("\nS6 global RI index:\n"); print(as.data.frame(s6 %>% select(contrast, n1, n2, delta_index, p, FDR)), row.names = FALSE)

## ── Figures ──────────────────────────────────────────────────────────────────
p_volc <- s1 %>% filter(!is.na(p)) %>%
  mutate(sig = ifelse(FDR < 0.05, "FDR < 0.05", ifelse(p < 0.05, "p < 0.05 (not FDR)", "n.s."))) %>%
  ggplot(aes(deltaPSI, -log10(p), colour = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  facet_wrap(~ contrast) +
  scale_colour_manual(values = c("FDR < 0.05" = "#C0392B", "p < 0.05 (not FDR)" = "#E59866", "n.s." = "grey70")) +
  labs(x = expression(Delta * "PSI (group 1 - group 2)"), y = expression(-log[10] * "(p)"), colour = NULL,
       title = "snoRNA-containing RI events: per-event Wilcoxon (libraries, lanes merged)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "S1_volcano_RI_events.png"), p_volc, width = 9, height = 7, dpi = 150)

p_idx <- ri_index %>% mutate(group = factor(group, levels = c("Healthy", "Normal", "Low", "High"))) %>%
  ggplot(aes(group, RI_index)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(x = NULL, y = "Global RI index (mean z-scored PSI)",
       title = "Global retention of snoRNA-containing introns per library") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "S6_global_RI_index.png"), p_idx, width = 6, height = 4.5, dpi = 150)

cat("\nDone. Tables in", OUT_DIR, "| figures in", FIG_DIR, "\n")

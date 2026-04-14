library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)

## 1) Load PSI (RI_events_with_snoRNA.tsv)
psi <- read_tsv("results/RI_events_with_snoRNA.tsv", show_col_types = FALSE)

## 2) Load SraRunTable.csv (metadata from GEO/SRA)
meta_raw <- read_csv("metadata/SraRunTable.csv", show_col_types = FALSE)

## 3) Build meta: sample + group
meta <- meta_raw %>%
  transmute(
    sample = Run,
    group = case_when(
      str_detect(tolower(condition), "healthy") ~ "Healthy",
      str_detect(tolower(condition), "normal")  ~ "Normal",
      str_detect(tolower(condition), "low")     ~ "Low",
      str_detect(tolower(condition), "high")    ~ "High",
      str_detect(tolower(grade), "low")         ~ "Low",
      str_detect(tolower(grade), "high")        ~ "High",
      str_detect(tolower(tissue), "adjacent") & str_detect(tolower(tissue), "normal") ~ "Normal",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct()

table(meta$group, useNA = "ifany")

## 4) Make long table (PSI + group)
psi_long <- psi %>%
  pivot_longer(
    cols = -c(eventID, snoRNAName),
    names_to = "sample",
    values_to = "PSI"
  ) %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(group))

n_distinct(psi_long$eventID)
table(psi_long$group)

## 5) Run differential for 4 contrasts
contrasts <- tibble(
  contrast = c("Low_vs_Healthy", "High_vs_Healthy", "Normal_vs_Healthy", "Low_vs_High"),
  g1 = c("Low", "High", "Normal", "Low"),
  g2 = c("Healthy", "Healthy", "Healthy", "High")
)

run_wilcox <- function(df, g1, g2) {
  df2 <- df %>% filter(group %in% c(g1, g2)) %>% drop_na(PSI)
  if (n_distinct(df2$group) < 2) return(tibble(p = NA_real_, deltaPSI = NA_real_))
  tibble(
    p = suppressWarnings(wilcox.test(PSI ~ group, data = df2)$p.value),
    deltaPSI = mean(df2$PSI[df2$group == g1], na.rm = TRUE) -
      mean(df2$PSI[df2$group == g2], na.rm = TRUE)
  )
}

event_results <- contrasts %>%
  pmap_dfr(function(contrast, g1, g2) {
    psi_long %>%
      group_by(eventID, snoRNAName) %>%
      group_modify(~ run_wilcox(.x, g1, g2)) %>%
      ungroup() %>%
      mutate(
        contrast = contrast,
        FDR = p.adjust(p, method = "BH")
      ) %>%
      relocate(contrast, eventID, snoRNAName, deltaPSI, p, FDR)
  })

write_tsv(event_results, "results/event_dPSI_results.tsv")

## 6) Quick checks
event_results %>%
  group_by(contrast) %>%
  summarise(total_events = n(), .groups = "drop")

event_results %>%
  mutate(sig = FDR < 0.05 & abs(deltaPSI) >= 0.10) %>%
  group_by(contrast) %>%
  summarise(
    sig_events = sum(sig, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  )

event_results %>%
  group_by(contrast) %>%
  summarise(
    p_median = median(p, na.rm = TRUE),
    p_min = min(p, na.rm = TRUE),
    .groups = "drop"
  )

event_results %>%
  filter(contrast == "High_vs_Healthy") %>%
  arrange(FDR) %>%
  select(eventID, deltaPSI, FDR) %>%
  head(10)

event_results %>%
  mutate(sig = FDR < 0.05 & abs(deltaPSI) >= 0.05) %>%
  group_by(contrast) %>%
  summarise(sig = sum(sig, na.rm = TRUE))

event_results %>%
  group_by(contrast) %>%
  summarise(
    n_FDR_lt_0.05 = sum(FDR < 0.05, na.rm = TRUE),
    n_FDR_lt_0.01 = sum(FDR < 0.01, na.rm = TRUE),
    .groups = "drop"
  )

snorna_summary <- event_results %>%
  mutate(sig = FDR < 0.05 & abs(deltaPSI) >= 0.10) %>%
  group_by(contrast, snoRNAName) %>%
  summarise(
    n_sig = sum(sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))

write_tsv(snorna_summary, "results/snorna_summary.tsv")

head(snorna_summary, 15)

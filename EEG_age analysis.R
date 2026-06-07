library(tidyverse)
library(rstatix)

base_dir <- "D:/EEG"

methods <- c("dSPM", "eLORETA")
conn_types <- c("imcoh", "wpli2_debiased")
bands <- c("theta", "alpha", "beta")

demo <- read_csv(file.path(base_dir, "Participants_LEMON.csv")) %>%
  rename(subject = `...1`) %>%
  mutate(
    subject = str_remove(subject, "^sub-"),
    subject = str_trim(subject),
    
    age_low = as.numeric(str_extract(Age, "^\\d+")),
    age_high = as.numeric(str_extract(Age, "\\d+$")),
    age_mid = (age_low + age_high) / 2,
    
    age_group_main = if_else(age_mid <= 60, "young", "old"),
    
    age_group_ablation = case_when(
      age_mid < 35 ~ "young",
      age_mid > 60 ~ "old",
      TRUE ~ NA_character_
    )
  )

#connectivity
read_conn <- function(inv_method) {
  path <- file.path(
    base_dir,
    paste0("connectivity_", inv_method),
    "stats",
    "bootstrap_connectivity_summary.tsv"
  )
  
  read_tsv(path) %>%
    mutate(
      inverse_method = inv_method,
      subject = str_remove(subject, "^sub-"),
      subject = str_trim(subject)
    )
}

conn <- bind_rows(read_conn("dSPM"), read_conn("eLORETA"))

conn2 <- conn %>% left_join(demo, by = "subject")

conn2 %>% count(age_group_main)

age_tests_main <- conn2 %>%
  filter(
    band %in% bands,
    method %in% conn_types,
    !is.na(age_group_main)
  ) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_main == "young"),
    n_old = sum(age_group_main == "old"),
    mean_young = mean(bootstrap_mean[age_group_main == "young"], na.rm = TRUE),
    mean_old = mean(bootstrap_mean[age_group_main == "old"], na.rm = TRUE),
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(bootstrap_mean ~ age_group_main, exact = FALSE)$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    comparison = "old_vs_young_main",
    p_fdr = p.adjust(p, method = "BH")
  ) %>%
  arrange(p_fdr)

age_tests_main

write_csv(age_tests_main, file.path(base_dir, "age_tests_main_FDR.csv"))

#inspect significant results
age_tests_main %>%
  select(inverse_method, method, band, task, p, p_fdr)

age_tests_main %>% filter(p_fdr < 0.05)

#ablation split
age_tests_ablation <- conn2 %>%
  filter(band %in% bands, method %in% conn_types, !is.na(age_group_ablation)) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_ablation == "young"),
    n_old = sum(age_group_ablation == "old"),
    
    mean_young = mean(
      bootstrap_mean[age_group_ablation == "young"], na.rm = TRUE),
    
    mean_old = mean(
      bootstrap_mean[age_group_ablation == "old"], na.rm = TRUE),
    
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(
        bootstrap_mean ~ age_group_ablation,
        exact = FALSE
      )$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    comparison = "old_vs_young_ablation",
    p_fdr = p.adjust(p, method = "BH")
  ) %>%
  arrange(p_fdr)

age_tests_ablation

age_tests_ablation %>% filter(p_fdr < 0.05)

arrange(age_tests_ablation, p)

ggplot(conn2 %>%
    filter(band == "alpha", method == "imcoh"),
  aes(age_group_ablation, bootstrap_mean, fill = age_group_ablation)) +
  geom_boxplot() + facet_grid(inverse_method ~ task) + theme_bw()

ggplot(conn2 %>%
    filter(band == "alpha", method == "imcoh"),
  aes(age_mid, bootstrap_mean)) + geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + facet_grid(inverse_method ~ task) + theme_bw()

#eyeopen vs eyeclosed
eyes_tests <- conn2 %>%
  filter(band %in% bands, method %in% conn_types) %>%
  select(subject, inverse_method, method, band, task, bootstrap_mean) %>%
  pivot_wider(names_from = task, values_from = bootstrap_mean) %>%
  filter(!is.na(EyesOpen), !is.na(EyesClosed)) %>%
  group_by(inverse_method, method, band) %>%
  summarise(
    p = wilcox.test(EyesOpen, EyesClosed, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop") %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  arrange(p_fdr)

eyes_tests

eyes_summary <- conn2 %>%
  filter(band == "alpha", method %in% c("imcoh", "wpli2_debiased")) %>%
  group_by(inverse_method, method, task) %>%
  summarise(
    mean_conn = mean(bootstrap_mean, na.rm = TRUE),
    sd_conn = sd(bootstrap_mean, na.rm = TRUE),
    .groups = "drop"
  )

eyes_summary

alpha_data <- conn2 %>%
  filter(band == "alpha", method == "wpli2_debiased", inverse_method == "dSPM")

ggplot(alpha_data, aes(task, bootstrap_mean, group = subject)) +
  geom_line(alpha = 0.3) + geom_point(size = 2) + theme_bw()

delta_alpha <- conn2 %>%
  filter(band == "alpha", method == "wpli2_debiased", inverse_method == "dSPM") %>%
  select(subject, age_mid, task, bootstrap_mean) %>%
  pivot_wider(names_from = task, values_from = bootstrap_mean) %>%
  mutate(delta = EyesClosed - EyesOpen)

ggplot(delta_alpha, aes(age_mid, delta)) + geom_point() +
  geom_smooth(method = "lm") + theme_bw() + labs(
    y = "Alpha Connectivity Difference\n(EyesClosed - EyesOpen)")

#power
read_power <- function(inv_method) {
  path <- file.path( base_dir, paste0("power_", inv_method), "stats",
                     "band_power_summary.tsv")
  
  read_tsv(path) %>%
    mutate(inverse_method = inv_method, subject = str_remove(subject, "^sub-"),
           subject = str_trim(subject))
}

power <- bind_rows(read_power("dSPM"), read_power("eLORETA"))

names(power)
head(power)

power2 <- power %>% left_join(demo, by = "subject")

power2 %>% count(age_group_main)

#age analysis for power
power_age_main <- power2 %>%
  filter( band %in% bands, !is.na(age_group_main)) %>%
  group_by(inverse_method, band, task) %>%
  summarise(
    n_young = sum(age_group_main == "young"),
    n_old = sum(age_group_main == "old"),
    
    mean_young = mean(global_log_power[age_group_main == "young"], na.rm = TRUE),
    
    mean_old = mean(global_log_power[age_group_main == "old"], na.rm = TRUE),
    
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(global_log_power ~ age_group_main, exact = FALSE)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(comparison = "power_old_vs_young_main", 
         p_fdr = p.adjust(p, method = "BH")) %>%
  arrange(p_fdr)

power_age_main

#ablation analysis for power
power_age_ablation <- power2 %>%
  filter(band %in% bands, !is.na(age_group_ablation)) %>%
  group_by(inverse_method, band, task) %>%
  summarise(
    n_young = sum(age_group_ablation == "young"),
    n_old = sum(age_group_ablation == "old"),
    
    mean_young = mean(global_log_power[age_group_ablation == "young"], na.rm = TRUE),
    
    mean_old = mean(global_log_power[age_group_ablation == "old"], na.rm = TRUE),
    
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(global_log_power ~ age_group_ablation, exact = FALSE)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(comparison = "power_old_vs_young_ablation", 
         p_fdr = p.adjust(p, method = "BH")) %>%
  arrange(p_fdr)

power_age_ablation

#eyeclosed vs eyeopen for power
power_eyes_tests <- power2 %>%
  filter(band %in% bands) %>%
  select(subject, inverse_method, band, task, global_log_power) %>%
  pivot_wider(names_from = task, values_from = global_log_power) %>%
  filter(!is.na(EyesOpen), !is.na(EyesClosed)) %>%
  group_by(inverse_method, band) %>%
  summarise(
    mean_EyesOpen = mean(EyesOpen, na.rm = TRUE),
    mean_EyesClosed = mean(EyesClosed, na.rm = TRUE),
    
    p = wilcox.test(EyesOpen, EyesClosed, paired = TRUE, exact = FALSE)$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(comparison = "power_EyesOpen_vs_EyesClosed", p_fdr = p.adjust(p, method = "BH")) %>%
  arrange(p_fdr)

power_eyes_tests

# stratified 
age_tests_main_stratified <- conn2 %>%
  filter(
    band %in% bands,
    method %in% conn_types,
    !is.na(age_group_main)
  ) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_main == "young"),
    n_old = sum(age_group_main == "old"),
    mean_young = mean(bootstrap_mean[age_group_main == "young"], na.rm = TRUE),
    mean_old = mean(bootstrap_mean[age_group_main == "old"], na.rm = TRUE),
    
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(bootstrap_mean ~ age_group_main, exact = FALSE)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  group_by(inverse_method, method) %>%
  mutate(p_fdr_stratified = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fdr_stratified)

age_tests_main_stratified

age_tests_ablation_stratified <- conn2 %>%
  filter(
    band %in% bands, method %in% conn_types, !is.na(age_group_ablation)) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_ablation == "young"),
    n_old = sum(age_group_ablation == "old"),
    mean_young = mean(bootstrap_mean[age_group_ablation == "young"], na.rm = TRUE),
    mean_old = mean(bootstrap_mean[age_group_ablation == "old"], na.rm = TRUE),
    
    p = ifelse(
      n_young > 0 & n_old > 0,
      wilcox.test(bootstrap_mean ~ age_group_ablation, exact = FALSE)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  group_by(inverse_method, method) %>%
  mutate(p_fdr_stratified = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fdr_stratified)

age_tests_ablation_stratified

# ks test
ks_age_main <- conn2 %>%
  filter(band %in% bands, method %in% conn_types, !is.na(age_group_main)) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_main == "young"),
    n_old = sum(age_group_main == "old"),
    ks_p = ifelse( n_young > 0 & n_old > 0, ks.test(
      bootstrap_mean[age_group_main == "young"],
      bootstrap_mean[age_group_main == "old"])$p.value, NA_real_),
    .groups = "drop") %>%
  group_by(inverse_method, method) %>%
  mutate(ks_p_fdr = p.adjust(ks_p, method = "BH")) %>%
  ungroup() %>%
  arrange(ks_p_fdr)

ks_age_main

ks_age_ablation <- conn2 %>%
  filter(band %in% bands, method %in% conn_types, !is.na(age_group_ablation)) %>%
  group_by(inverse_method, method, band, task) %>%
  summarise(
    n_young = sum(age_group_ablation == "young"),
    n_old = sum(age_group_ablation == "old"),
    ks_p = ifelse(
      n_young > 0 & n_old > 0,
      ks.test(
        bootstrap_mean[age_group_ablation == "young"],
        bootstrap_mean[age_group_ablation == "old"]
      )$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  group_by(inverse_method, method) %>%
  mutate(ks_p_fdr = p.adjust(ks_p, method = "BH")) %>%
  ungroup() %>%
  arrange(ks_p_fdr)

ks_age_ablation

ggplot(
  conn2 %>%
    filter(band == "alpha", method == "imcoh", !is.na(age_group_ablation)),
  aes(x = age_group_ablation, y = bootstrap_mean)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.7) +
  facet_grid(inverse_method ~ task) + theme_bw() +
  labs(x = "Age Group", y = "Connectivity",
    title = "Alpha imcoh Connectivity: Young <35 vs Old >60"
  )

anova_data <- conn2 %>%
  filter(band %in% bands, method %in% conn_types, !is.na(age_group_main))

anova_model <- aov(
  bootstrap_mean ~ age_group_main * band * task * inverse_method * method,
  data = anova_data
)

summary(anova_model)

ancova_model <- aov(
  bootstrap_mean ~ age_mid * band * task * inverse_method * method,
  data = anova_data
)

summary(ancova_model)

# dSPM vs eLORETA using paired testing Wilcoxon-Signed rank test
dspm_vs_eloreta <- conn2 %>%
  filter(band %in% bands, method %in% conn_types) %>%
  select(subject, task, band, method, inverse_method, bootstrap_mean) %>%
  pivot_wider(names_from = inverse_method, values_from = bootstrap_mean) %>%
  filter(!is.na(dSPM), !is.na(eLORETA)) %>%
  group_by(method, band, task) %>%
  summarise(mean_dSPM = mean(dSPM, na.rm = TRUE),
            mean_eLORETA = mean(eLORETA, na.rm = TRUE),
            p = wilcox.test(dSPM, eLORETA, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(method) %>%
  mutate(p_fdr_stratified = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fdr_stratified)

dspm_vs_eloreta

# WPLI-debiased vs IMCOH using paired testing Wilcoxon-Signed rank test
wpli_vs_imcoh <- conn2 %>%
  filter(band %in% bands, method %in% conn_types) %>%
  select(subject, task, band, inverse_method, method, bootstrap_mean) %>%
  pivot_wider(names_from = method, values_from = bootstrap_mean) %>%
  filter(!is.na(wpli2_debiased), !is.na(imcoh)) %>%
  group_by(inverse_method, band, task) %>%
  summarise(mean_wpli = mean(wpli2_debiased, na.rm = TRUE),
            mean_imcoh = mean(imcoh, na.rm = TRUE),
            p = wilcox.test(wpli2_debiased, imcoh, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(inverse_method) %>%
  mutate(p_fdr_stratified = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fdr_stratified)

wpli_vs_imcoh

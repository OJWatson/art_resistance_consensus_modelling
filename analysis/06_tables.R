# Get data for generating tables
ls <- grep("csv", list.files("analysis/data/raw/sim_outputs", recursive = TRUE, full.names = TRUE), value = TRUE)
dat <- lapply(ls, read.csv, stringsAsFactors = FALSE)
for(x in seq_along(ls)) {
  dat[[x]][["rep"]] <- NULL
  dat[[x]][["X"]] <- NULL
  dat[[x]]$MODEL <- lapply(stringr::str_split(ls, "/"),"[[", 5)[[x]]
  dat[[x]]$SCENARIO <- substr(lapply(stringr::str_split(ls, "/"),"[[", 6)[[x]], 1, 2)
}

df <- dplyr::bind_rows(dat)
df$MODEL <- factor(c("MORU", "PSU", "Imperial")[match(df$MODEL,c("moru", "psu", "imperial"))], levels=c("MORU", "PSU", "Imperial"))
df_cache <- df[df$STARTING_PARTNER_DRUG_FREQ < 1, ]

# Downstream table generation packages
library(tidyverse)
library(broom)
library(survival)

## summary functions
survival_replace <- function(x) {

  to_replace <- sum(x == 40)
  if(to_replace == 0) {
    return(x)
  } else if(to_replace == length(x)) {
    return(rep(999, length(x)))
  } else {

  df <- data.frame("ty" = x)
  df$completed <- as.integer(df$ty < 40)
  df_surv <- with(df, survival::Surv(ty, event = completed))
  weib.fit <- survival::survreg(df_surv~1, df, dist="weibull")
  x[x==40] <- tail(predict(weib.fit,data.frame("df_surv" = 1),p = 0:99/100, type = "quantile"), to_replace)
  return(x)
  }
}

perc_diff_func <- function(x, y) {

  wd <- DescTools::MeanDiffCI(x, y, method = "basic")
  mean <- paste0(sprintf("%.1f", round(100*wd[1] / mean(x),1)), "%")
  ci <- paste(paste0(sprintf("%.1f", round(100*wd[2:3] / mean(x),1), "%")), collapse = ", ")
  return(paste0(mean, " [", ci, "]"))

}

x_to_z_if_y <- function(x, y = 45, z = "40+ ") {if(x == y) {z} else {sprintf("%.1f", x)}}

# set all values of NA to 40. These will then be imputed using Weibull survival
df  <- df_cache
df[,grep("TIME", names(df))][which(df[,grep("TIME", names(df))] >= 40, arr.ind = TRUE)] <- NA
df[,grep("TIME", names(df))][which(is.na(df[,grep("TIME", names(df))]), arr.ind = TRUE)] <- 40

# Group 580Y milestones
test <- df %>% group_by(MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ, PFPR, TREATMENT_COVERAGE) %>%
  mutate(TIME_TO_1p_580Y = survival_replace(TIME_TO_1p_580Y),
         TIME_TO_10p_580Y = survival_replace(TIME_TO_10p_580Y),
         TIME_TO_25p_580Y = survival_replace(TIME_TO_25p_580Y))

# first table summary of times with 0 PD
tbl_a <- test %>%
  filter(TREATMENT_COVERAGE == 0.4 & STARTING_PARTNER_DRUG_FREQ == 0) %>%
  group_by(MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ) %>%
  summarise(TY_25 = paste0(x_to_z_if_y(median(TIME_TO_25p_580Y)),
                           " [",
                           paste(x_to_z_if_y(quantile(TIME_TO_25p_580Y, 0.25)),
                                 x_to_z_if_y(quantile(TIME_TO_25p_580Y, 0.75)),
                                 sep = ", "),
                           "]")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A4", "DHA-PPQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A5", "ASAQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A6", "AL")) %>%
  arrange(desc(SCENARIO)) %>%
  tidyr::pivot_wider(names_from = STARTING_PARTNER_DRUG_FREQ, values_from = TY_25,
                     names_prefix = "PD_")

# Split out 0 PD into new column
test_formed <- left_join(
  test %>% filter(TREATMENT_COVERAGE == 0.4 & STARTING_PARTNER_DRUG_FREQ != 0.0) %>% mutate(rep = 1:100),
  test %>% filter(TREATMENT_COVERAGE == 0.4 & STARTING_PARTNER_DRUG_FREQ == 0.0) %>% mutate(rep = 1:100) %>%
    rename(TIME_TO_25p_580Y_0 = TIME_TO_25p_580Y) %>% ungroup %>%
    select(PFPR, TREATMENT_COVERAGE, SCENARIO, MODEL, TIME_TO_25p_580Y_0, rep)
  )

# Create second table for percentage difference
tbl_b <- test_formed %>%
  filter(TREATMENT_COVERAGE == 0.4) %>%
  group_by(MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ) %>%
  summarise(TY_25_DIFF = perc_diff_func(TIME_TO_25p_580Y_0, TIME_TO_25p_580Y)) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A4", "DHA-PPQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A5", "ASAQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A6", "AL")) %>%
  arrange(desc(SCENARIO)) %>%
  tidyr::pivot_wider(names_from = STARTING_PARTNER_DRUG_FREQ, values_from = TY_25_DIFF,
                     names_prefix = "PD_DIFF")

# Lastly get the gradients
grads <-   test %>%
  filter(TREATMENT_COVERAGE == 0.4) %>%
  ungroup %>%
  select(MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ, TIME_TO_25p_580Y) %>%
  nest(data = -(1:2)) %>%
  mutate(
    fit = map(data, ~ lm(TIME_TO_25p_580Y ~ STARTING_PARTNER_DRUG_FREQ, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term != "(Intercept)") %>%
  select(MODEL, SCENARIO, estimate) %>%
  mutate(estimate = sprintf("%.1f", round(estimate/-10, 2))) %>%
  rename(`Years Lost` = estimate) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A4", "DHA-PPQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A5", "ASAQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A6", "AL"))

# Table 1 Overall
to_save <- left_join(left_join(tbl_a, tbl_b), grads)
write.table(to_save, file = "analysis/tables/tbl1.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# Supp Table
# first table summary of times with 0 PD
supp_table <- test %>%
  filter(TREATMENT_COVERAGE == 0.4 & STARTING_PARTNER_DRUG_FREQ == 0) %>%
  group_by(MODEL, PFPR, SCENARIO, STARTING_PARTNER_DRUG_FREQ) %>%
  summarise(TY_25 = paste0(x_to_z_if_y(median(TIME_TO_25p_580Y)),
                           " [",
                           paste(x_to_z_if_y(quantile(TIME_TO_25p_580Y, 0.25)),
                                 x_to_z_if_y(quantile(TIME_TO_25p_580Y, 0.75)),
                                 sep = ", "),
                           "]")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A4", "DHA-PPQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A5", "ASAQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A6", "AL")) %>%
  arrange(desc(SCENARIO)) %>%
  tidyr::pivot_wider(names_from = c("STARTING_PARTNER_DRUG_FREQ","PFPR"), values_from = TY_25,
                     names_prefix = "PD_")
names(supp_table)[3:6] <- paste0("PfPR = ", c(1, 5, 10, 20), "%")
write.table(supp_table, file = "analysis/tables/supp_table_pfpr.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# years lost per change
yls <- test %>%
  filter(TREATMENT_COVERAGE == 0.4) %>%
  group_by(MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ) %>%
  summarise(TY_25 = median(TIME_TO_25p_580Y)) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A4", "DHA-PPQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A5", "ASAQ")) %>%
  mutate(SCENARIO = replace(SCENARIO, SCENARIO=="A6", "AL")) %>%
  arrange(desc(SCENARIO)) %>%
  tidyr::pivot_wider(names_from = STARTING_PARTNER_DRUG_FREQ, values_from = TY_25,
                     names_prefix = "PD_")
yls$PD_0.01 <- yls$PD_0 - yls$PD_0.01
yls$PD_0.1 <- yls$PD_0 - yls$PD_0.1
yls$PD_0.25 <- yls$PD_0 - yls$PD_0.25
yls$PD_0.5 <- yls$PD_0 - yls$PD_0.5
yls$subsequent <- (yls$PD_0.5 - yls$PD_0.1)/4

# supp figure for PfPR
library(scales)

# modify x axis data for easy log plotting
modify_dat <- function(x){
  x$STARTING_PARTNER_DRUG_FREQ[x$STARTING_PARTNER_DRUG_FREQ==0] <- 0.001
  return(x)
}

supp_plot_pfpr <- test %>%
  mutate(SCENARIO = replace(SCENARIO, seq_along(SCENARIO), c("DHA-PPQ", "ASAQ", "AL")[match(SCENARIO, c("A4","A5","A6"))])) %>%
  mutate(PFPR = factor(paste0(PFPR, "%"), levels = paste0(c(1,5,10,20), "%"))) %>%
  filter(TREATMENT_COVERAGE == 0.4) %>% modify_dat %>%
  ggplot(aes(STARTING_PARTNER_DRUG_FREQ, TIME_TO_25p_580Y, color = as.factor(PFPR))) +
  geom_jitter(alpha  = 0.5, width = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(MODEL~SCENARIO, scales = "free") +
  ggpubr::theme_pubclean() + theme(axis.line = element_line()) +
  scale_color_viridis_d(name = "PfPR:") +
  scale_x_log10(breaks = c(0.001, (c(0.01,0.1,0.25,0.5))),labels = c("0", (c(0.01,0.1,0.25,0.5)))) +
  ylab("Years to 0.25 580Y") +
  xlab("Starting Partner Drug Resistance Frequency") +
  theme(legend.key = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "grey"),
        panel.border = element_rect(colour = "grey", fill = NA))
save_figs("supp_pfpr", supp_plot_pfpr, width = 9, height = 6)

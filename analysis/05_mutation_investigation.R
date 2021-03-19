library(tidyverse)

## MUTATION ANALYSIS BY IMPERIAL

# Bring in raw simulations for Imperial Model
df <- readRDS("analysis/data/raw/mutation_checks/imperial_mutation_analysis.rds")

# relabel loci
df <- df %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L1", "pfcrt76T")) %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L2", "pfmdr186Y")) %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L3", "pfmdr1184F")) %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L4", "pfmdr1CNV")) %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L5", "pfk13580Y")) %>%
  mutate(Drug_Locus = replace(Drug_Locus, Drug_Locus == "L6", "pfpm2CNV")) %>%
  mutate(Drug_Locus = factor(Drug_Locus, levels = c("pfcrt76T", "pfmdr186Y","pfmdr1184F","pfmdr1CNV","pfk13580Y","pfpm2CNV"))) %>%
  mutate(selection = as.integer(Drug_Locus)>4) %>%
  mutate(Time = replace(Time, Time == "TIME_TO_1p", "Years to 0.01 Frequency")) %>%
  mutate(Time = replace(Time, Time == "TIME_TO_10p", "Years to 0.1 Frequency")) %>%
  mutate(Time = replace(Time, Time == "TIME_TO_25p", "Years to 0.25 Frequency")) %>%
  mutate(Mu = paste("Mu = ", Mu))

# create plot
gg_mu <- ggplot(df, aes(x = Drug_Locus, y = value, color = selection)) +
  facet_grid(Mu~Time) +
  geom_boxplot() +
  ylab("Years") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Locus Under Selection from DHA-PPQ:", option = "B", end = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.line = element_line(), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white"), legend.key = element_rect(fill = "white"),
        axis.title.x = element_blank())

save_figs("mutation_investigation", gg_mu, width = 10, height = 6)

## COTRANSMISSION ANALYSIS BY IMPERIAL

# Bring in raw simulations for Imperial Model
df <- read.csv("analysis/data/raw/sa_cotransmission/SA_COTRANS_IMPERIAL_magenta_1.3.0_20210121.csv")
df$SPZ_label <- c("1 sporozoite per infection", "2 sporozoites per infection")[match(df$SPZ, c(0.01, 0.2))]
df$FIRST_LINE <- factor(c("DHA-PPQ", "ASAQ", "AL")[df$FIRST_LINE+1], levels = c("DHA-PPQ", "ASAQ", "AL"))

# create plot
gg_cotrans <- df %>% ggplot(aes(x = TIME_TO_1p_580Y, y = as.factor(STARTING_PARTNER_DRUG_FREQ), color = SPZ_label)) +
  geom_boxplot(notch = FALSE) + facet_wrap(~FIRST_LINE, scales = "free_x") +
  #ylab("Starting Partner Drug Resistance Frequency") +
  ylab("Starting Partner Drug \nResistance Frequency") +
  xlab("Years till 0.01 580Y Frequency") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Cotransmission:", option = "B", end = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.line = element_line(), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white"), legend.key = element_rect(fill = "white"))

gg_cotrans_10 <- df %>% ggplot(aes(x = TIME_TO_10p_580Y, y = as.factor(STARTING_PARTNER_DRUG_FREQ), color = SPZ_label)) +
  geom_boxplot(notch = FALSE) + facet_wrap(~FIRST_LINE, scales = "free_x") +
  ylab("Starting Partner Drug \nResistance Frequency") +
  xlab("Years till 0.1 580Y Frequency") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Cotransmission:", option = "B", end = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.line = element_line(), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white"), legend.key = element_rect(fill = "white"))

gg_cotrans_25 <- df %>% ggplot(aes(x = TIME_TO_25p_580Y, y = as.factor(STARTING_PARTNER_DRUG_FREQ), color = SPZ_label)) +
  geom_boxplot(notch = FALSE) + facet_wrap(~FIRST_LINE, scales = "free_x") +
  ##ylab("Starting Partner Drug Resistance Frequency") +
  ylab("Starting Partner Drug \nResistance Frequency") +
  xlab("Years till 0.25 580Y Frequency") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Cotransmission:", option = "B", end = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.line = element_line(), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white"), legend.key = element_rect(fill = "white"),
        axis.title.y = element_blank())


gg_cotrans <- cowplot::plot_grid(gg_cotrans, gg_cotrans_10, gg_cotrans_25, ncol = 1)

save_figs("cotransmission_investigation", gg_cotrans, width = 10, height = 10)

## Load in the drug by genotype table
library(tidyverse)
tbls <- "analysis/data/uncertainty_drug_tables/EfficacyTable_0.csv"
tbls <- vapply(0:9, function(x) {gsub("0", x, tbls)}, character(1))

tbls <- lapply(tbls, readr::read_csv, locale = readr::locale(encoding = "UTF-16"))
for(i in 1:10) {tbls[[i]]$tbl <- i}

tbls <- do.call(rbind, tbls)
gg1 <- tbls %>% select(ID, Genotype, ASAQ, AL, DHAPPQ, tbl) %>%
  group_by(ID, Genotype) %>%
  summarise(across(ASAQ:DHAPPQ, .fns = list(
    min = ~quantile(.x, 0.025),
    max = ~quantile(.x, 0.975),
    med = ~quantile(.x, 0.5)
  )
  )) %>%
  pivot_longer(ASAQ_min:DHAPPQ_med, names_to = c("drug", "stat"), names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  ggplot(aes(ID, med, ymin = min, ymax = max)) +
  geom_point(size = 0.5) +
  geom_errorbar(size = 0.5) +
  facet_wrap(~drug) +
  ylab("28-Day Treatment Efficacy") +
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1),
        #axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none")
save_figs("B1_DxG_tables_style1", gg1, width = 12, height = 6)

gg2 <- tbls %>%
  mutate(
    Genotype = gsub("--", "0", Genotype),
    Genotype = gsub("x", "", Genotype),
    Genotype = gsub("X", "", Genotype),
    Genotype = gsub("NFNF", "NF1", Genotype),
    Genotype = gsub("YYYY", "YY1", Genotype),
    Genotype = gsub("YFYF", "YF1", Genotype),
    Genotype = gsub("NYNY", "NY1", Genotype)
  ) %>%
  select(Genotype, ASAQ, AL, DHAPPQ, tbl) %>%
  unique %>%
  pivot_longer(ASAQ:DHAPPQ, names_to = c("drug"), names_pattern = "(.*)") %>%
  arrange(drug) %>%
  mutate(ID = rep(1:64, n()/64)) %>%
  ggplot(aes(ID, value, group = tbl, color = as.character(tbl))) +
  geom_point(size = 0.5) +
  facet_wrap(~drug) +
  ylab("28-Day Treatment Efficacy") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Drug By Genotype Efficacy Table:") +
  theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1),
        #axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white")) +
  xlab("Genotype ID")
save_figs("B1_DxG_tables_style2", gg2, width = 12, height = 6)

### and the B2 analysis table

tbls <- "analysis/data/uncertainty_drug_tables/EfficacyTable_equal.csv"
tbls <- readr::read_csv(tbls, locale = readr::locale(encoding = "UTF-16"))

gg3 <- tbls %>%
  mutate(tbl = 1) %>%
  mutate(
    Genotype = gsub("--", "0", Genotype),
    Genotype = gsub("x", "", Genotype),
    Genotype = gsub("X", "", Genotype),
    Genotype = gsub("NFNF", "NF1", Genotype),
    Genotype = gsub("YYYY", "YY1", Genotype),
    Genotype = gsub("YFYF", "YF1", Genotype),
    Genotype = gsub("NYNY", "NY1", Genotype)
  ) %>%
  select(Genotype, ASAQ, AL, DHAPPQ, tbl) %>%
  unique %>%
  pivot_longer(ASAQ:DHAPPQ, names_to = c("drug"), names_pattern = "(.*)") %>%
  arrange(drug) %>%
  mutate(ID = rep(1:64, n()/64)) %>%
  ggplot(aes(ID, value, group = drug, color = as.character(drug))) +
  geom_point(size = 2) +
  facet_wrap(~drug) +
  ylab("28-Day Treatment Efficacy") +
  ggpubr::theme_pubclean() +
  scale_color_viridis_d(name = "Drug By Genotype Efficacy Table:") +
  theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1),
        #axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white")) +
  xlab("Genotype ID")
save_figs("B2_DxG_table", gg3, width = 12, height = 6)


# Plot our default efficacy table:
what <- magenta::drug_table
what$Genotype <- drug_table$Genotype

gens <- c("KNY0C0",
          "KNY0Y0",
          "KNY0C1", "TNY0C0","KNF1C0",
          "KNY0Y1",  "TYY0Y0", "KNF1Y0")
al_gens <- gens[c(1,2,5,8)]
dhappq_gens <- gens[c(1,2,3,6)]
asaq_gens <- gens[c(1,2,4,7)]
what <- what %>%
  filter(Genotype %in% gens) %>%
  select(DHAPPQ, ASAQ, AL, Genotype)
what <- what[match(gens, what$Genotype), ]


what$type <- "Wildtype"
what$type[2] <- "Wildtype + 580Y"
what$type[3:5] <- "Starting PD Mutant"
what$type[6:8] <- "Starting PD Mutant + 580Y"
what$type <- factor(what$type, levels = c("Wildtype", "Wildtype + 580Y", "Starting PD Mutant", "Starting PD Mutant + 580Y"))
what[c(4,5,7,8),1] <- NA
what[c(3,5,6,8),2] <- NA
what[c(3,4,6,7),3] <- NA
gg4 <- what %>% pivot_longer(DHAPPQ:AL) %>%
  mutate(Genotype = gsub("1$", "2", Genotype),
         Genotype = gsub("0$", "1", Genotype)) %>%
  ggplot(aes(type, value, color = name, label = Genotype )) +
  ggrepel::geom_label_repel(show.legend = FALSE) +
  geom_jitter(width = 0.02) +
  xlab("") +
  theme_bw() +
  ylab("28 Day Treatment Efficacy") +
  scale_color_discrete(name = "Drug")
save_figs("deafult_efficacy_table", gg4, 8, 6)

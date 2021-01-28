library(tidyverse)
library(ggridges)

# ------------------------------------------------------------------------------
## 1. Get the data -------------------------------------------------------------
# ------------------------------------------------------------------------------

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
df <- df[df$STARTING_PARTNER_DRUG_FREQ < 1, ]

# set all times that werent reached to grrreater than 40 for presentation ease
df[,grep("TIME", names(df))][which(df[,grep("TIME", names(df))] >= 40, arr.ind = TRUE)] <- NA
df[,grep("TIME", names(df))][which(is.na(df[,grep("TIME", names(df))]), arr.ind = TRUE)] <- 45

# ------------------------------------------------------------------------------
## 2. Summaries of -------------------------------------------------------------
# ------------------------------------------------------------------------------

# check plot
summaries <- lapply(c("A4", NA,"A5", NA, "A6"), function(x) {

  if(is.na(x)) {
    return(NULL)
  }

  title <- c("Dihydroartemisinin-Piperaquine", "Artesunate-Amodiaquine", "Artemether-Lumefantrine")[match(x, c("A4", "A5", "A6"))]
  title <- paste("Scenarios with", title, "used as the first-line therapy")
  gg <- df[df$SCENARIO==x, ] %>%
    rowwise() %>%
  mutate(PFPR = if(nchar(as.character(PFPR))==1){ paste0("0",as.character(PFPR)) } else { as.character(PFPR) }) %>%
    mutate(PFPR = paste0("PfPR = ", PFPR, "%")) %>%
    mutate(TREATMENT_COVERAGE = scales::percent(TREATMENT_COVERAGE, prefix = " Coverage = ", suffix = "%")) %>%
    ungroup %>%
    group_by(STARTING_PARTNER_DRUG_FREQ, MODEL, PFPR, TREATMENT_COVERAGE) %>%
    summarise(y0 = min(TIME_TO_25p_580Y),
  y25 = quantile(TIME_TO_25p_580Y, 0.25),
  y50 = median(TIME_TO_25p_580Y),
  y75 = quantile(TIME_TO_25p_580Y, 0.75),
  y100 = max(TIME_TO_25p_580Y))
  gg[gg == 45] <- 40


  gg_plot <- ggplot(gg, aes(x=as.factor(STARTING_PARTNER_DRUG_FREQ),
                                   fill=as.factor(MODEL),
                                   group = interaction(STARTING_PARTNER_DRUG_FREQ, MODEL))) +
  facet_grid(PFPR~TREATMENT_COVERAGE) +
  geom_hline(yintercept = 10, linetype="dotted") +
  theme_bw() +
    xlab("Starting partner drug resistance frequency") +
    ylab("Years for 580Y frequency to reach 0.25") +
  geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100,),
               width=0.4, position = position_dodge(width=0.7), stat="identity") +
  scale_fill_manual(values = c("#7f9fad","#dda7c7","#ffd27f"), name = "") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=10),
          legend.position = "top",
          legend.key.size = unit(18, "pt"),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(title) + ylim(c(0,40))

  return(gg_plot)
})

# individual plot saves
leg <- cowplot::get_legend(summaries[[1]] + theme(legend.position = "right"))
save_figs("dha_ppq", cowplot::plot_grid(summaries[[1]] + theme(legend.position = "none"), leg, rel_widths = c(1,0.1)), width = 12, height = 6)
save_figs("asaq", cowplot::plot_grid(summaries[[3]] + theme(legend.position = "none"), leg, rel_widths = c(1,0.1)), width = 12, height = 6)
save_figs("al", cowplot::plot_grid(summaries[[5]] + theme(legend.position = "none"), leg, rel_widths = c(1,0.1)), width = 12, height = 6)

# alltogether
summaries[[1]] <- summaries[[1]] + theme(legend.position = "none")
summaries[[3]] <- summaries[[3]] + theme(axis.title.y = element_blank(), legend.position = "none")
summaries[[5]] <- summaries[[5]] + theme(legend.position = "none", axis.title.y = element_blank())

leg <- cowplot::get_legend(summaries[[1]] + theme(legend.position = "top", legend.key.size = unit(36, "pt")))
gg <- cowplot::plot_grid(plotlist = summaries, ncol = 5, rel_widths = c(1,0.05,1,0.05,1))
gg2 <- cowplot::plot_grid(leg, gg, ncol =1, rel_heights = c(0.1,1))
save_figs("all_drugs", gg2, width = 20, height = 10)

# ------------------------------------------------------------------------------
## 3. Median Times -------------------------------------------------------------
# ------------------------------------------------------------------------------

# set NAs to 45 so that we get the correct median
med_df <- group_by(df, MODEL, SCENARIO, STARTING_PARTNER_DRUG_FREQ, PFPR, TREATMENT_COVERAGE) %>%
  summarise_all(.funs = median, na.rm=TRUE) %>%
  ungroup()

# and change them to NA for plotting
med_df[med_df==45] <- NA
rel_25p <- lapply(c("A4", NA, "A5", NA, "A6"), function(x) {

  if(is.na(x)) {
    return(NULL)
  }

  title <- c("Dihydroartemisinin-Piperaquine", "Artesunate-Amodiaquine", "Artemether-Lumefantrine")[match(x, c("A4", "A5", "A6"))]

  gg_df <- med_df %>%
    filter(SCENARIO == x) %>%
    rowwise() %>%
    mutate(PFPR = if(nchar(as.character(PFPR))==1){ paste0("0",as.character(PFPR)) } else { as.character(PFPR) }) %>%
    mutate(PFPR = paste0("PfPR = ", PFPR, "%")) %>%
    mutate(TREATMENT_COVERAGE = scales::percent(TREATMENT_COVERAGE, prefix = "Coverage = ", suffix = "%"))
  gg_df$TIME_TO_25p_580Y_TEXT <- vapply(gg_df$TIME_TO_25p_580Y, function(x) {
    ret <- as.character(round(x))
    if(is.na(ret)) {
      return("40+")
    } else {
      return(ret)
    }
  }, character(1))

  gg <- gg_df %>%
    ggplot(aes(x=MODEL, y=as.factor(STARTING_PARTNER_DRUG_FREQ), fill=TIME_TO_25p_580Y)) +
    geom_tile(color="white") +
    geom_text(aes(label = TIME_TO_25p_580Y_TEXT), color = "white") +
    ylab("Starting partner drug resistance frequency") +
    facet_grid(PFPR~TREATMENT_COVERAGE) +
    scale_fill_viridis_c(option = "B",
                         limits = c(0, 40.1),
                         direction = -1,
                         end = 0.85,
                        guide = "colourbar") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    labs(fill = "Years for 580Y frequency to reach 0.25") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          panel.grid.minor = element_line(colour="grey", size=0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.box="horizontal",
          strip.background = element_blank(),
          strip.text = element_text(size=10),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 35)) +
    ggtitle(title)

  return(gg)
})

leg <- cowplot::get_legend(rel_25p[[1]])
rel_25p[[1]] <- rel_25p[[1]] + theme(legend.position = "none")
rel_25p[[3]] <- rel_25p[[3]] + theme(axis.title.y = element_blank(), legend.position = "none")
rel_25p[[5]] <- rel_25p[[5]] + theme(legend.position = "none",axis.title.y = element_blank())

gg <- cowplot::plot_grid(plotlist = rel_25p, ncol = 5,
                         rel_widths = c(1,0.05,1,0.05,1),
                         labels = c("a", "", "b", "", "c"))
gg2 <- cowplot::plot_grid(leg, gg, ncol =1, rel_heights = c(0.1,1))
save_figs(name = "median_times_tiled", gg2, height = 9, width = 15, root = "analysis/plots")


# ------------------------------------------------------------------------------
## 3. Selctions, s -------------------------------------------------------------
# ------------------------------------------------------------------------------

df <- do.call(rbind, dat)
df$MODEL[df$MODEL=="PSU_CIDD"] <- "PSU"
df$MODEL <- factor(df$MODEL, levels=c("MORU", "PSU", "IMPERIAL"))
df <- df[df$STARTING_PARTNER_DRUG_FREQ < 1, ]
df$TIME_TO_1p_580Y[df$TIME_TO_1p_580Y>40] <- NA
df$TIME_TO_10p_580Y[df$TIME_TO_10p_580Y>40] <- NA
df$TIME_TO_25p_580Y[df$TIME_TO_25p_580Y>40] <- NA

df$s01025 <- (log(0.25/(1-0.25)) - log(0.10/(1-0.1))) / ((df$TIME_TO_25p_580Y-df$TIME_TO_10p_580Y)*6)
df$s00101 <- (log(0.1/(1-0.1)) - log(0.01/(1-0.01))) / ((df$TIME_TO_10p_580Y-df$TIME_TO_1p_580Y)*6)

ridge_25p <- lapply(c("A4", NA, "A5", NA, "A6"), function(x) {

  if(is.na(x)) {
    return(NULL)
  }

  gg <- df[df$SCENARIO == x, ] %>%
    filter(STARTING_PARTNER_DRUG_FREQ != 0) %>%
    filter(SCENARIO == x) %>%
    rowwise() %>%
    mutate(PFPR = if(nchar(as.character(PFPR))==1){ paste0("0",as.character(PFPR)) } else { as.character(PFPR) }) %>%
    mutate(PFPR = paste0("PfPR = ", PFPR, "%")) %>%
    mutate(TREATMENT_COVERAGE = scales::percent(TREATMENT_COVERAGE, prefix = "Coverage = ", suffix = "%"))

ggplot(gg, aes(x=s01025, y = as.factor(STARTING_PARTNER_DRUG_FREQ),
                     fill=as.factor(MODEL))) +
  geom_density_ridges(alpha=0.8, scale = 1) +
  facet_grid(PFPR~TREATMENT_COVERAGE, scales = "free_x") +
  theme_bw() +
  ylab("Starting partner drug resistance (%)") +
  xlab("Selection Coefficient from 0.1 to 0.25 580Y frequency") +
  scale_fill_manual(values = c("#7f9fad","#dda7c7","#ffd27f"), name = "") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10),
        legend.position = "top",
        legend.key.size = unit(36, "pt")) +
  ggtitle(x)

})

leg <- cowplot::get_legend(ridge_25p[[1]])
ridge_25p[[1]] <- ridge_25p[[1]] + theme(legend.position = "none")
ridge_25p[[3]] <- ridge_25p[[3]] + theme(axis.title.y = element_blank(), legend.position = "none")
ridge_25p[[5]] <- ridge_25p[[5]] + theme(legend.position = "none",axis.title.y = element_blank())

gg <- cowplot::plot_grid(plotlist = ridge_25p, ncol = 5, rel_widths = c(1,0.05,1,0.05,1))
gg2 <- cowplot::plot_grid(leg, gg, ncol =1, rel_heights = c(0.1,1))
save_figs("selection_ridges", gg2, height = 10, width = 20)


# ------------------------------------------------------------------------------
## 4. Selctions, l -------------------------------------------------------------
# ------------------------------------------------------------------------------

s1_gg <- df %>% group_by(SCENARIO, MODEL, PFPR, TREATMENT_COVERAGE, STARTING_PARTNER_DRUG_FREQ ) %>%
  summarise(s = median(s00101, na.rm = TRUE)) %>%
  mutate(drug = factor(c("DHA-PPQ", "ASAQ", "AL")[match(SCENARIO, c("A4", "A5", "A6"))], levels = c("DHA-PPQ", "ASAQ", "AL"))) %>%
  mutate(pfpr_print = factor(paste("PfPR:", PFPR, "%"), levels = unique(paste("PfPR:", PFPR, "%")))) %>%
  mutate(TREATMENT_COVERAGE = scales::percent(TREATMENT_COVERAGE, prefix = "Coverage = ", suffix = "%")) %>%
  filter(TREATMENT_COVERAGE != "Coverage = 40%") %>%
  ggplot(aes(x = STARTING_PARTNER_DRUG_FREQ, y = s, color = MODEL, shape = as.factor(TREATMENT_COVERAGE))) +
  geom_point() +
  facet_grid(drug~pfpr_print) +
  geom_line(aes(linetype = as.factor(TREATMENT_COVERAGE)), lwd = 1) +
  scale_x_sqrt(breaks = unique(df$STARTING_PARTNER_DRUG_FREQ), labels = unique(df$STARTING_PARTNER_DRUG_FREQ)) +
  scale_y_log10(breaks = c(0, 0.02,0.05,0.1,0.2), labels = c(0.0, 0.02, 0.05, 0.1, 0.2)) +
  geom_hline(yintercept = c(0, 0.02,0.05,0.1), lwd = 0) +
  scale_color_manual(values = c("#577988","#c870a4","#ffb733"), name = "Model") +
  scale_linetype(name = "Treatment Coverage") +
  scale_shape(name = "Treatment Coverage") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10)) +
  xlab("Starting Partner Drug Resistance Frequency\n") +
  ylab("\nSelection Coefficient, s (0.01 0.1)")


s2_gg <- df %>% group_by(SCENARIO, MODEL, PFPR, TREATMENT_COVERAGE, STARTING_PARTNER_DRUG_FREQ ) %>%
  summarise(s = median(s01025, na.rm = TRUE)) %>%
  mutate(drug = factor(c("DHA-PPQ", "ASAQ", "AL")[match(SCENARIO, c("A4", "A5", "A6"))], levels = c("DHA-PPQ", "ASAQ", "AL"))) %>%
  mutate(pfpr_print = factor(paste("PfPR:", PFPR, "%"), levels = unique(paste("PfPR:", PFPR, "%")))) %>%
  mutate(TREATMENT_COVERAGE = scales::percent(TREATMENT_COVERAGE, prefix = "Coverage = ", suffix = "%")) %>%
  filter(TREATMENT_COVERAGE != "Coverage = 40%") %>%
  ggplot(aes(x = STARTING_PARTNER_DRUG_FREQ, y = s, color = MODEL, shape = as.factor(TREATMENT_COVERAGE))) +
  geom_point() +
  facet_grid(drug~pfpr_print) +
  geom_line(aes(linetype = as.factor(TREATMENT_COVERAGE)), lwd = 1) +
  scale_x_sqrt(breaks = unique(df$STARTING_PARTNER_DRUG_FREQ), labels = unique(df$STARTING_PARTNER_DRUG_FREQ)) +
  scale_y_log10(breaks = c(0, 0.02,0.05,0.1,0.2), labels = c(0.0, 0.02, 0.05, 0.1, 0.2)) +
  geom_hline(yintercept = c(0, 0.02,0.05,0.1), lwd = 0) +
  scale_color_manual(values = c("#577988","#c870a4","#ffb733"), name = "Model") +
  scale_linetype(name = "Treatment Coverage") +
  scale_shape(name = "Treatment Coverage") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10)) +
  xlab("Starting Partner Drug Resistance Frequency\n") +
  ylab("\nSelection Coefficient, s (0.1 0.25)")

save_figs("selection_coeffs", cowplot::plot_grid(s1_gg, s2_gg, ncol = 1, labels = "auto"), width = 10, height = 10)

# ------------------------------------------------------------------------------
## 4. Selctions, h -------------------------------------------------------------
# ------------------------------------------------------------------------------

s3_gg <- df %>%
  mutate(drug = factor(c("DHA-PPQ", "ASAQ", "AL")[match(SCENARIO, c("A4", "A5", "A6"))], levels = c("DHA-PPQ", "ASAQ", "AL"))) %>%
  ggpubr::ggdensity(x = "s01025", fill = "MODEL", color = "MODEL",
                    add.params = list(size = 1.25, linetype = "dashed", show.legend = FALSE),
                    add = "median", y = "..scaled..", facet.by = "drug", alpha = 0.2) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(),
        strip.background = element_rect(fill = "white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(color = "grey", size = 0.25),
        panel.grid.minor = element_blank(),
        ) +
  scale_fill_manual(values = c("#577988","#c870a4","#e69500"), name = "Model:") +
  scale_color_manual(values = c("#577988","#c870a4","#e69500"), name = "Model:") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,0.3,0.05)) +
  geom_hline(yintercept = 1.1, lwd = 0) +
  theme(panel.background = element_rect(colour = "grey", linetype = "solid"))+
  xlab(expression(Selection~Coefficient~(s[0.1-0.25]))) +
  ylab("Scaled Density") +
  guides(color = guide_legend(override.aes = list(colour = NULL, lwd = 0)),
         fill = guide_legend(override.aes = list(alpha = 1)))

save_figs("selection_densities", s3_gg, width = 12, height = 5)

# ------------------------------------------------------------------------------
## 4. Alignments ---------------------------------------------------------------
# ------------------------------------------------------------------------------

alignment_moru <- read_csv('analysis/data/raw/figure1/moru/A4_MORU_rev17.csv') %>% mutate(MODEL = "MORU")
alignment_psu <- read_csv('analysis/data/raw/figure1/psu/A4_PSU_mu_0p001983_20200229.csv') %>% mutate(MODEL = "PSU")
alignment_imperial = read_csv('analysis/data/raw/figure1/imperial/A4_IMPERIAL_magenta_1.3.0_2020022.csv')[,-1] %>% mutate(MODEL = "Imperial")
alignment_df <- dplyr::bind_rows(aligment_moru, alignment_psu, alignment_imperial) %>%
  mutate(MODEL = factor(MODEL, levels=c("MORU", "PSU", "Imperial")))

alignment_gg <- alignment_df %>%
  filter(TREATMENT_COVERAGE == 0.4 & PFPR == 10 & STARTING_PARTNER_DRUG_FREQ == 0) %>%
  ggpubr::ggdensity(x = "TIME_TO_1p_580Y", fill = "MODEL", color = "MODEL",
                    add.params = list(size = 1.25, linetype = "dashed", show.legend = FALSE),
                    add = "median", y = "..scaled..",alpha = 0.2) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(),
        strip.background = element_rect(fill = "white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(color = "grey", size = 0.25),
        panel.grid.minor = element_blank(),
  ) +
  scale_fill_manual(values = c("#577988","#c870a4","#e69500"), name = "Model:") +
  scale_color_manual(values = c("#577988","#c870a4","#e69500"), name = "Model:") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1,13,2)) +
  geom_hline(yintercept = 1.1, lwd = 0) +
  theme(panel.background = element_rect(colour = "grey", linetype = "solid"))+
  xlab("Time to 0.01 580Y") +
  ylab("Scaled Density") +
  guides(color = guide_legend(override.aes = list(colour = NULL, lwd = 0)),
         fill = guide_legend(override.aes = list(alpha = 1)))

save_figs("alignment", alignment_gg, width = 6, height = 3)

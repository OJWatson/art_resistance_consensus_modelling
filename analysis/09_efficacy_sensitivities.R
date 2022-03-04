library(tidyverse)
devtools::load_all()

# ------------------------------------------------------------------------------
## 1. Get the data -------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----------------
# B1
# ----------------

ls <- grep("csv", list.files("analysis/data/raw/B1_results/", recursive = TRUE, full.names = TRUE), value = TRUE)
ls <- grep("0p10", ls, value = TRUE)

dat <- lapply(ls, read.csv, stringsAsFactors = FALSE)
for(x in seq_along(ls)) {
  dat[[x]][["rep"]] <- NULL
  dat[[x]][["X"]] <- NULL
  dat[[x]]$MODEL <- lapply(stringr::str_split(ls[x], "/"),"[[", 6)[[1]]
  dat[[x]]$MODEL <- factor(c("MORU", "PSU", "Imperial")[match(dat[[x]]$MODEL,c("MORU", "PSU", "IMPERIAL"))], levels=c("MORU", "PSU", "Imperial"))
  tbl_ns <- grep("Table", stringr::str_split(lapply(stringr::str_split(ls[x], "/"),"[[", 7)[[1]], "_")[[1]])+1
  dat[[x]]$TABLE <- as.numeric(stringr::str_split(lapply(stringr::str_split(ls[x], "/"),"[[", 7)[[1]], "_")[[1]][tbl_ns])
}

df <- dplyr::bind_rows(dat)
df <- df[df$STARTING_PARTNER_DRUG_FREQ < 1, ] # PSU Had a an extra starting partner drug scenario here

# set all times that werent reached to grrreater than 40 for presentation ease
df[,grep("TIME", names(df))][which(df[,grep("TIME", names(df))] >= 40, arr.ind = TRUE)] <- NA
df[,grep("TIME", names(df))][which(is.na(df[,grep("TIME", names(df))]), arr.ind = TRUE)] <- 45

# set up drugs
df$FIRST_LINE <- factor(c("DHA-PPQ", "ASAQ", "AL")[df$FIRST_LINE+1], levels = c("DHA-PPQ", "ASAQ", "AL"))
b1 <- df

# ----------------
# B2
# ----------------

ls <- grep("csv", list.files("analysis/data/raw/B2_results/", recursive = TRUE, full.names = TRUE), value = TRUE)

dat <- lapply(ls, read.csv, stringsAsFactors = FALSE)
for(x in seq_along(ls)) {
  dat[[x]][["rep"]] <- NULL
  dat[[x]][["X"]] <- NULL
  dat[[x]]$MODEL <- lapply(stringr::str_split(ls[x], "/"),"[[", 6)[[1]]
  dat[[x]]$MODEL <- factor(c("MORU", "PSU", "Imperial")[match(dat[[x]]$MODEL,c("MORU", "PSU", "IMPERIAL"))], levels=c("MORU", "PSU", "Imperial"))
  dat[[x]]$TABLE <- ((x)-1) %% 2
}

df <- dplyr::bind_rows(dat)
df <- df[df$STARTING_PARTNER_DRUG_FREQ < 1, ] # PSU Had a an extra starting partner drug scenario here

# set all times that werent reached to grrreater than 40 for presentation ease
df[,grep("TIME", names(df))][which(df[,grep("TIME", names(df))] >= 40, arr.ind = TRUE)] <- NA
df[,grep("TIME", names(df))][which(is.na(df[,grep("TIME", names(df))]), arr.ind = TRUE)] <- 45

# set up drugs
df$FIRST_LINE <- factor(c("DHA-PPQ", "ASAQ", "AL")[df$FIRST_LINE+1], levels = c("DHA-PPQ", "ASAQ", "AL"))
b2 <- df


# ------------------------------------------------------------------------------
## 2. Summary plots -------------------------------------------------------------
# ------------------------------------------------------------------------------

# B1

# check plot
summaries_b1 <- function(df, models = c("MORU", "PSU", "Imperial")) {

  plots <- lapply(models, function(x) {

    gg <- df %>% filter(MODEL == x) %>%
      group_by(STARTING_PARTNER_DRUG_FREQ, MODEL, FIRST_LINE, TABLE) %>%
      summarise(y0 = min(TIME_TO_25p_580Y),
                y25 = quantile(TIME_TO_25p_580Y, 0.25),
                y50 = median(TIME_TO_25p_580Y),
                y75 = quantile(TIME_TO_25p_580Y, 0.75),
                y100 = max(TIME_TO_25p_580Y))
    suppressWarnings(gg[gg > 40] <- 40)

    gg %>%
      mutate(TABLE = interaction(as.factor(TABLE), FIRST_LINE, STARTING_PARTNER_DRUG_FREQ)) %>%
      mutate(STARTING_PARTNER_DRUG_FREQ = paste0("Starting Partner Drug \nResistance Frequency = ", STARTING_PARTNER_DRUG_FREQ)) %>%
      ggplot(aes(x=fct_reorder(TABLE, desc(y50)),
                 fill=as.factor(FIRST_LINE))) +
      facet_grid(MODEL~STARTING_PARTNER_DRUG_FREQ, scales = "free_x") +
      geom_hline(yintercept = 10, linetype="dotted") +
      theme_bw() +
      xlab("Genotype by Drug Efficacy Table Ranked by Median Selection Times") +
      ylab("Years for 580Y frequency to reach 0.25") +
      # geom_errorbar(aes(ymin = y0, y = y50, ymax = y100),
      #              width=0.4, position = position_dodge(width=0.7), stat="identity") +
      geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                   #color = NA,
                   width=0.4, position = position_dodge(width=0.7), stat="identity") +
      # geom_point(aes(y = y50),shape = "-",position = position_dodge(width=0.7), stat="identity", size = 3) +
      scale_fill_manual(values = c("#7f9fad","#dda7c7","#ffd27f"), name = "") +
      theme(strip.background = element_blank(),
            strip.text = element_text(size=12),
            legend.position = "top",
            legend.key.size = unit(18, "pt"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank()) +
      ylim(c(0,40)) +
      theme(legend.key.width =  unit(30, "pt"), legend.text = element_text(size = 12),
            legend.position = "top", text = element_text(size = 12), strip.text.x = element_text(size = 10))

  })




  return(plots)
}

b1_summaries <- summaries_b1(b1)
b1_plot <- cowplot::plot_grid(
  cowplot::get_legend(b1_summaries[[2]]),
  b1_summaries[[1]] + xlab("") + ylab("") + theme(axis.ticks.x = element_blank(), legend.position = "none", axis.title.x = element_blank()),
  b1_summaries[[2]] + xlab("") + theme(axis.ticks.x = element_blank(), legend.position = "none", axis.title.x = element_blank(), strip.text.x = element_blank()),
  b1_summaries[[3]] + xlab("") + ylab("") + theme(axis.ticks.x = element_blank(), legend.position = "none", strip.text.x = element_blank()),
  ncol = 1,
  rel_heights = c(0.1, 1.05,1.05, 1))
save_figs("B1", b1_plot, width = 12, height = 10)



# B2
summaries_b2 <- function(df, tbls = 0, models = c("MORU", "PSU", "Imperial")) {

  plots <- lapply(models, function(x) {

    gg <- df %>% filter(MODEL == x) %>% filter(TABLE %in% tbls) %>%
      group_by(STARTING_PARTNER_DRUG_FREQ, MODEL, FIRST_LINE, TABLE) %>%
      summarise(y0 = min(TIME_TO_25p_580Y),
                y25 = quantile(TIME_TO_25p_580Y, 0.25),
                y50 = median(TIME_TO_25p_580Y),
                y75 = quantile(TIME_TO_25p_580Y, 0.75),
                y100 = max(TIME_TO_25p_580Y)) %>%
      mutate(tablelabel = interaction(paste0("v",as.factor(TABLE+1)), FIRST_LINE, STARTING_PARTNER_DRUG_FREQ))

    suppressWarnings(gg[gg > 40] <- 40)

    gg %>%
      mutate(STARTING_PARTNER_DRUG_FREQ = paste0("Starting Partner Drug \nResistance Frequency = ", STARTING_PARTNER_DRUG_FREQ)) %>%
      ggplot(aes(x=fct_reorder(tablelabel, desc(y50)),
                 fill=as.factor(FIRST_LINE))) +
      facet_grid(MODEL~STARTING_PARTNER_DRUG_FREQ, scales = "free_x") +
      geom_hline(yintercept = 10, linetype="dotted") +
      theme_bw() +
      xlab("Genotype by Drug Efficacy Table Ranked by Median Selection Times") +
      ylab("Years for 580Y frequency to reach 0.25") +
      # geom_errorbar(aes(ymin = y0, y = y50, ymax = y100),
      #              width=0.4, position = position_dodge(width=0.7), stat="identity") +
      geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                   #color = NA,
                   width=0.4, position = position_dodge(width=0.7), stat="identity") +
      # geom_point(aes(y = y50),shape = "-",position = position_dodge(width=0.7), stat="identity", size = 3) +
      scale_fill_manual(values = c("#7f9fad","#dda7c7","#ffd27f"), name = "") +
      theme(strip.background = element_blank(),
            strip.text = element_text(size=12),
            legend.position = "top",
            legend.key.size = unit(18, "pt"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ylim(c(0,40)) +
      theme(legend.key.width =  unit(30, "pt"), legend.text = element_text(size = 12),
            legend.position = "top", text = element_text(size = 12), strip.text.x = element_text(size = 10))

  })




  return(plots)
}


b2_summaries <- summaries_b2(b2, tbls = 0)
b2_plot <- cowplot::plot_grid(
  cowplot::get_legend(b2_summaries[[2]]),
  b2_summaries[[1]] + xlab("") + ylab("\n") + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()),
  b2_summaries[[2]] + xlab("") + ylab("Years for 580Y frequency \nto reach 0.25") +
    theme( legend.position = "none", axis.title.x = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()),
  b2_summaries[[3]]  + xlab("") + ylab("\n") + theme( legend.position = "none", strip.text.x = element_blank(), axis.text.x = element_blank()),
  ncol = 1,
  rel_heights = c(0.1, 0.95, 0.8, 0.8))
b2_plot
save_figs("B2", b2_plot, width = 12, height = 7)

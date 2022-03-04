library(tidyverse)
devtools::load_all()

# ------------------------------------------------------------------------------
## 1. Get the data -------------------------------------------------------------
# ------------------------------------------------------------------------------

# Get the drug lookups for PSU data
drug_table <- magenta::drug_table
genotypeInfo <- read.table('scripts/genotypeID.txt', header = TRUE);
genotypeInfo$Shortname <- gsub("--", "0", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("x", "", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("X", "", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("NFNF", "NF1", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("YYYY", "YY1", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("YFYF", "YF1", genotypeInfo$Shortname)
genotypeInfo$Shortname <- gsub("NYNY", "NY1", genotypeInfo$Shortname)
genotypeInfo$ID <- genotypeInfo$ID + 1
genotypeInfo$imp_id <- match(genotypeInfo$Shortname, drug_table$Genotype)

# ----------------
# A4-A6 PSU LINEAGES
# ----------------

ls <- grep("PSU.*monthly", list.files("analysis/data/raw/sim_lineages/", recursive = TRUE, full.names = TRUE), value = TRUE)

get_psu_lineages <- function(x, rep) {

  dat <- read.table(x, stringsAsFactors = FALSE)
  drugs <- gsub("(.*)(_monthly_data)","\\1",unlist(strsplit(x, "/") %>% lapply("[[", 7)))

  monthly_data <- dat
  genotype_distribution <- monthly_data[,20:147]
  rowSum <- rowSums(genotype_distribution)
  for(i in seq_along(nrow(monthly_data))) {
    genotype_distribution[i,] <- genotype_distribution[i,]/rowSum[i]
  }

  # final 482 rows are the last 40 years from Nguyen
  gap <- diff(seq(0, 40, length.out = 482)+20)[1]
  names(genotype_distribution) <- genotypeInfo$Shortname
  genotype_distribution$Time <- seq(gap, nrow(genotype_distribution)*gap, gap)
  genotype_distribution$Time <- genotype_distribution$Time - max(genotype_distribution$Time) + 40
  res <- genotype_distribution %>% pivot_longer(cols = genotypeInfo$Shortname) %>%
    rename(Genotype = name, Strains = value)


  res2 <- res %>% group_by(Time,  Genotype) %>%
    summarise(Strains = sum(Strains)) %>%
    mutate(model = "PSU") %>%
    mutate(rep = rep) %>%
    mutate(drug = drugs)

  return(res2)

}

dat <- vector("list", length(ls))
for(i in seq_along(dat)) {
  dat[[i]] <- get_psu_lineages(ls[i], i %% 100)
}

psu_dat <- do.call(rbind, dat)

psu_dat_x <- lapply(c("DHAPPQ", "ASAQ", "AL"), function(x) {
  y <- psu_dat  %>% filter(drug == x)
  y %>% filter(Genotype %in% unique(y$Genotype[which(y$Strains>0.05)]))
})

pd_ar_plot <- function(x, pd, drug, model = "PSU"){

  x %>%
    mutate(art = grepl("Y.$", Genotype)) %>%
    mutate(pd = grepl(pd, Genotype)) %>%
    group_by(rep, Time) %>%
    summarise(pd = sum(Strains[pd]),
              art = sum(Strains[art]),
              wt = 1 - pd - art) %>%
    ggplot(aes(Time, pd - art, group = rep))  +
    geom_line(alpha = 0.2) +
    geom_line(aes(Time, med),
              data = x %>%
                mutate(art = grepl("Y.$", Genotype)) %>%
                mutate(pd = grepl(pd, Genotype)) %>%
                group_by(rep, Time) %>%
                summarise(pd = sum(Strains[pd]),
                          art = sum(Strains[art]),
                          wt = 1 - pd - art) %>%
                group_by(Time) %>%
                summarise(med = median(pd - art)), color = "red", lwd = 2,
              inherit.aes = FALSE) +
    ylab("580Y Frequency - Maximally Resistant \nPartner Drug Frequency") +
    theme_bw() +
    ggtitle(paste0(model, ": ", drug)) +
    xlim(c(-1, 40))

}

psu_gg_1 <- pd_ar_plot(psu_dat %>% filter(drug == "DHAPPQ"), "2$", "DHA-PPQ") + ylab("580Y Frequency - Plasmepsin Gene \nAmplifcation Frequency")
psu_gg_2 <- pd_ar_plot(psu_dat %>% filter(drug == "ASAQ"), "TYY|KYY", "ASAQ") + ylab("580Y Frequency - pfmdr1  86Y and Y184 \nHaplotype Frequency")
psu_gg_3 <- pd_ar_plot(psu_dat %>% filter(drug == "AL"), "KNF|TNF", "AL") + ylab("580Y Frequency - pfmdr1  N86 and 184F \nHaplotype Frequency")

psu_lineage_gg <- cowplot::plot_grid(psu_gg_1, psu_gg_2, psu_gg_3, ncol = 1)


# -----------------------------------------------------------------------------

dhappq_res <- read.csv("analysis/data/raw/sim_lineages/IMPERIAL/DHAPPQ_monthly.csv")
asaq_res <- read.csv("analysis/data/raw/sim_lineages/IMPERIAL/ASAQ_monthly.csv")
al_res <- read.csv("analysis/data/raw/sim_lineages/IMPERIAL/AL_monthly.csv")

dhappq_res <- dhappq_res %>%
  mutate(drug = "DHAPPQ") %>%
  mutate(Genotype = gsub("1$", "2", Genotype),
         Genotype = gsub("0$", "1", Genotype)) %>%
  mutate(Time = Time - 20)
asaq_res <- asaq_res %>%
  mutate(drug = "ASAQ") %>%
  mutate(Genotype = gsub("1$", "2", Genotype),
         Genotype = gsub("0$", "1", Genotype)) %>%
  mutate(Time = Time - 20)
al_res <- al_res %>%
  mutate(drug = "AL") %>%
  mutate(Genotype = gsub("1$", "2", Genotype),
         Genotype = gsub("0$", "1", Genotype)) %>%
  mutate(Time = Time - 20)

imperial_gg_1 <- pd_ar_plot(dhappq_res, "2$", "DHA-PPQ", "Imperial") + ylab("580Y Frequency - Plasmepsin Gene \nAmplifcation Frequency")
imperial_gg_2 <- pd_ar_plot(asaq_res, "TYY|KYY", "ASAQ", "Imperial") + ylab("580Y Frequency - pfmdr1  86Y and Y184 \nHaplotype Frequency")
imperial_gg_3 <- pd_ar_plot(al_res, "KNF|TNF", "AL", "Imperial") + ylab("580Y Frequency - pfmdr1  N86 and 184F \nHaplotype Frequency")

imperial_lineage_gg <- cowplot::plot_grid(imperial_gg_1, imperial_gg_2, imperial_gg_3, ncol = 1)

lineage_full_gg <- cowplot::plot_grid(psu_lineage_gg, NULL, imperial_lineage_gg, ncol = 3, rel_widths = c(1,0.1,1))
save_figs("580Y_vs_PD_lineages", lineage_full_gg, width = 12, height = 12)

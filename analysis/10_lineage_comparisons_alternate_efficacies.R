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
# B2 PSU LINEAGES
# ----------------

ls <- grep("monthly", list.files("analysis/data/raw/sim_lineages/", recursive = TRUE, full.names = TRUE), value = TRUE)

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
  message(i)
  dat[[i]] <- get_psu_lineages(ls[i], i %% 100)
}

psu_dat <- do.call(rbind, dat)

gens <- as.character(unique(psu_dat$Genotype[psu_dat$Strains>0.05]))
color <- sample(pals::stepped(length(gens)))
color[gens == "KNY0C1"] <- "black"
names(color) <- gens

psu_dat2 <- lapply(c("DHAPPQ", "ASAQ", "AL"), function(x) {
  y <- psu_dat %>% filter(drug == x)
  y <- y %>% filter(Genotype %in% unique(y$Genotype[which(y$Strains>0.05)])) %>%
    group_by(Time, Genotype, drug) %>%
    summarise(y = median(Strains),
              ymin = quantile(Strains, 0.025),
              ymax = quantile(Strains, 0.975))
})

psu_dat_x <- lapply(c("DHAPPQ", "ASAQ", "AL"), function(x) {
  y <- psu_dat  %>% filter(drug == x)
  y %>% filter(Genotype %in% unique(y$Genotype[which(y$Strains>0.05)]))
})


gg1 <- ggplot(psu_dat_x[[1]], aes(Time, Strains, color = Genotype)) +
  geom_line(aes(group = interaction(rep, Genotype, drug)), alpha = 0.2) +
  geom_line(aes(y = y, x = Time, color = Genotype), data = psu_dat2[[1]], inherit.aes = FALSE, lwd = 2) +
  theme_bw() +
  scale_color_manual(values = color[names(color) %in% psu_dat_x[[1]]$Genotype], drop = TRUE) +
  ggtitle("PSU: DHA-PPQ") +
  xlim(c(-5,40))
gg1

gg2 <- ggplot(psu_dat_x[[2]], aes(Time, Strains, color = Genotype)) +
  geom_line(aes(group = interaction(rep, Genotype, drug)), alpha = 0.2) +
  geom_line(aes(y = y, x = Time, color = Genotype), data = psu_dat2[[2]], inherit.aes = FALSE, lwd = 2) +
  theme_bw() +
  scale_color_manual(values = color[names(color) %in% psu_dat_x[[2]]$Genotype]) +
  ggtitle("PSU: ASAQ") +
  xlim(c(-5,40))

gg3 <- ggplot(psu_dat_x[[3]], aes(Time, Strains, color = Genotype)) +
  geom_line(aes(group = interaction(rep, Genotype, drug)), alpha = 0.2) +
  geom_line(aes(y = y, x = Time, color = Genotype), data = psu_dat2[[3]], inherit.aes = FALSE, lwd = 2) +
  theme_bw() +
  scale_color_manual(values = color[names(color) %in% psu_dat_x[[3]]$Genotype]) +
  ggtitle("PSU: AL") +
  xlim(c(-5,40))

psu_lineage_gg <- cowplot::plot_grid(gg1, gg2, gg3, ncol = 1)


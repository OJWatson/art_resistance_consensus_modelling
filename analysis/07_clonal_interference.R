## POST
make_lineage_plot <- function(model, drug) {

  df <- read.csv(paste0("analysis/data/raw/lineages/",model,"_",drug,".csv"))
  drug_name <- c("DHA-PPQ","AL","ASAQ")[match(drug, c("dhappq","al","asaq"))]

  fits <- unique(cbind(df$Lineage, df$Prob_LPF))
  fits <- fits[order(fits[,2]),]

  df$Lineage <- factor(df$Lineage, levels = fits[,1])

  df %>%
    filter(Strains != 0) %>%
    ggplot(aes(Time, Strains, fill = Prob_LPF, group = Lineage)) +
    geom_bar(position = "fill",stat = "identity", lwd = 0) +
    ggpubr::theme_pubclean() +
    scale_fill_gradient2(midpoint = min(fitness_dhappq$Prob_LPF) + (max(fitness_dhappq$Prob_LPF)-min(fitness_dhappq$Prob_LPF))/2,
                         name = paste0("Probaility of 28 day treatment failure with ", drug_name)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("\nTime (years)") + ylab("Strain Proportions\n") +
    theme(axis.line = element_line(), panel.grid = element_blank()) +
    scale_x_continuous(expand = c(0,0), limits = c(0,40)) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20))

}

psu_dhappq <- make_lineage_plot("psu","dhappq")
psu_asaq <- make_lineage_plot("psu","asaq")
psu_al <- make_lineage_plot("psu","al")
imperial_dhappq <- make_lineage_plot("imperial","dhappq")
imperial_asaq <- make_lineage_plot("imperial","asaq")
imperial_al <- make_lineage_plot("imperial","al")

psu <- cowplot::plot_grid(psu_dhappq + xlab(""), NA, psu_asaq+ xlab(""), NA, psu_al+ xlab(""),
                          ncol = 5, rel_widths = c(1,0.1,1,0.1,1))

imperial <- cowplot::plot_grid(imperial_dhappq+theme(legend.position = "none"), NA,
                               imperial_asaq+theme(legend.position = "none"), NA,
                               imperial_al+theme(legend.position = "none"),
                               ncol = 5, rel_widths = c(1,0.1,1,0.1,1))

clonal_interf <- cowplot::plot_grid(psu, imperial, ncol = 1, rel_heights = c(1,0.9), labels = "auto")
save_figs("clonal_interference", clonal_interf, height = 10, width = 16)

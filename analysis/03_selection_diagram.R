# deomstration of stages of drug resistance

ex <- read.csv("analysis/data/raw/figure1/imperial/A4_init_res_freq_0p25_tc_0p4_pfpr_PFPR05_C580Y_Imperial.csv")
ex$time <- seq(0, 40, length.out = nrow(ex))
df <- tidyr::pivot_longer(ex, cols = head(everything(),-1)) %>%
  mutate(rep = as.numeric(unlist(lapply(strsplit(name, ".", fixed = TRUE), function(x) {x[2]}))))
df$rep[is.na(df$rep)] <- 0

# make plot
gg <- ggplot(df %>% filter(rep < 10), aes(time, value, group = rep)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(time, value), data = df %>% group_by(time) %>% summarise(value = median(value)), inherit.aes = FALSE, lwd = 2) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line()) +
  ylab("Frequency of Resistance Mutations") +
  xlab("Years") +
  geom_segment(aes(x=x,xend=xend,y=y,yend=yend),
               data = data.frame(x = 0, xend = 5, y = 0.1, yend = 0.1),
               arrow = arrow(ends = "both",type = "closed", angle = 20, length = unit(10, "pt")), inherit.aes = FALSE) +
  geom_segment(aes(x=x,xend=xend,y=y,yend=yend),
               data = data.frame(x = 5, xend = 18, y = 0.5, yend = 0.5),
               arrow = arrow(ends = "both",type = "closed", angle = 20, length = unit(10, "pt")), inherit.aes = FALSE) +
  geom_segment(aes(x=x,xend=xend,y=y,yend=yend),
               data = data.frame(x = 18, xend = 30, y = 0.9, yend = 0.9),
               arrow = arrow(ends = "both",type = "closed", angle = 20, length = unit(10, "pt")), inherit.aes = FALSE) +
  geom_text(aes(x=x,y=y,label=label), data = data.frame(x = 2.5, y = 0.15,label="Emergence"),size=4,inherit.aes = FALSE) +
  geom_text(aes(x=x,y=y,label=label), data = data.frame(x = 11.5, y = 0.55,label="Establishment"),size=4,inherit.aes = FALSE) +
  geom_text(aes(x=x,y=y,label=label), data = data.frame(x = 24, y = 0.95,label="Fixation"),size=4,inherit.aes = FALSE)

save_figs("resistance_stages", gg, width = 8, height = 6)

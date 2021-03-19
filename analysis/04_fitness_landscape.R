## Load in the drug by genotype table
library(tidyverse)
drug_table <- read.csv("analysis/data/raw/drug_table.csv")

# Sort it into format to make a network
dt_df <- data.table::rbindlist(lapply(lapply(strsplit(drug_table$Genotype, ""),t),as.data.frame))
names(dt_df) <- c("pfcrt76T", "pfmdr186Y", "pfmdr1184F", "pfmdr1CNV", "pfk13580Y", "pfpm2CNV")
dt_df$pfcrt76T <- as.numeric(factor(dt_df$pfcrt76T)) - 1
dt_df$pfmdr186Y <- as.numeric(factor(dt_df$pfmdr186Y)) - 1
dt_df$pfmdr1184F <- as.numeric(factor(dt_df$pfmdr1184F, levels = c("Y", "F"))) - 1
dt_df$pfmdr1CNV <- as.numeric(factor(dt_df$pfmdr1CNV)) - 1
dt_df$pfk13580Y <- as.numeric(factor(dt_df$pfk13580Y)) - 1
dt_df$pfpm2CNV <- as.numeric(factor(dt_df$pfpm2CNV)) - 1

# Build the graph for this
binaries <- apply(dt_df[,1:6], 1, paste0, collapse = "")
monos <- unique(t(apply(which(as.matrix(dist(dt_df[,1:6],diag = TRUE)) == 1, arr.ind = TRUE), 1, sort)))
g <- igraph::make_graph(as.numeric(t(monos)))
g <- igraph::set.vertex.attribute(g, name = "genotype", value = paste0(substr(drug_table$Genotype,1,3),"\n",substr(drug_table$Genotype,4,6)))
g <- igraph::set.vertex.attribute(g, name = "DP_fitness", value = drug_table$DHAPPQ)
g <- igraph::set.vertex.attribute(g, name = "ASAQ_fitness", value = drug_table$ASAQ)
g <- igraph::set.vertex.attribute(g, name = "AL_fitness", value = drug_table$AL)

# Set up a color pallete to show it
library(RColorBrewer)
colfunc <- colorRampPalette(colors = brewer.pal(11, "rainbow")[3:9])
brs_dhappq <- seq(min(1-drug_table$DHAPPQ), max(1-drug_table$DHAPPQ), length.out = 50)
brs_al <- seq(min(1-drug_table$AL), max(1-drug_table$AL), length.out = 50)
brs_asaq <- seq(min(1-drug_table$DHAPPQ), max(1-drug_table$ASAQ), length.out = 50)
cols <- gplots::rich.colors(length(brs_dhappq))

# Create the DHA_PPQ graph
svg("analysis/plots/dha_ppq_network.svg", width = 18,height = 8)

plot(g, layout = igraph::layout.reingold.tilford(g),
     vertex.label = igraph::vertex_attr(g, "genotype"),
     vertex.color = cols[unlist(lapply(1-drug_table$DHAPPQ, function(x) {which.min(abs(x - brs_dhappq))}))],
     vertex.frame.color = "black",
     vertex.frame.cex = c(1,2)[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     edge.width = 1,
     edge.arrow.width = 0.1,
     vertex.size = 6.5,
     edge.arrow.size = 0.2,
     vertex.label.color = "#ffffff",
     vertex.size2 = 4,
     vertex.label.family = "Lucida Console",
     vertex.label.cex = 0.95,
     vertex.shape = c("square","circle")[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     asp = 0.345,
     margin = -0)

fields::image.plot(legend.only=T, zlim=range(brs_dhappq), col=cols,
                   legend.lab = "\n\nProbability of 28-Day DHA-PPQ Treatment Failure",
                   legend.mar = 3, horizontal = TRUE)
dev.off()

# Create the AL graph
svg("analysis/plots/al_network.svg", width = 18,height = 8)

plot(g, layout = igraph::layout.reingold.tilford(g),
     vertex.label = igraph::vertex_attr(g, "genotype"),
     vertex.color = cols[unlist(lapply(1-drug_table$AL, function(x) {which.min(abs(x - brs_dhappq))}))],
     vertex.frame.color = "black",
     vertex.frame.cex = c(1,2)[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     edge.width = 1,
     edge.arrow.width = 0.1,
     vertex.size = 6.5,
     edge.arrow.size = 0.2,
     vertex.label.color = "#ffffff",
     vertex.size2 = 4,
     vertex.label.cex = 0.9,
     vertex.shape = c("square","circle")[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     asp = 0.345,
     margin = -0)

fields::image.plot(legend.only=T, zlim=range(brs_dhappq), col=cols,
                   legend.lab = "\n\nProbability of 28-Day AL Treatment Failure",
                   legend.mar = 3, horizontal = TRUE)
dev.off()


## Alternate style

dt_df$dist_dhappq <- as.matrix(dist(dt_df[,1:6],diag = TRUE, method = "manhattan"))[,1]
dt_df$dist_asaq <- as.matrix(dist(dt_df[,1:6],diag = TRUE, method = "manhattan"))[,1]
dt_df$dist_al <- as.matrix(dist(dt_df[,1:6],diag = TRUE, method = "manhattan"))[,1]
dt_df$fitness_dhappq <-  as.numeric(drug_table$DHAPPQ)
dt_df$fitness_asaq <-  as.numeric(drug_table$ASAQ)
dt_df$fitness_al <-  as.numeric(drug_table$AL)

g1 <- ggplot(dt_df, aes(dist_dhappq, fitness_dhappq)) + geom_point() + ylim(c(0.4,1)) +
        ylab("Treatment Efficacy by DHA-PPQ") +
        xlab("Distance from Wild Type") +
        ggpubr::theme_pubclean() +
        theme(axis.line = element_line())

g2 <- ggplot(dt_df, aes(dist_asaq, fitness_asaq)) + geom_point() + ylim(c(0.4,1)) +
        ylab("Treatment Efficacy by ASAQ") +
        xlab("Distance from Wild Type") +
        ggpubr::theme_pubclean() +
        theme(axis.line = element_line())

g3 <- ggplot(dt_df, aes(dist_al, fitness_al)) + geom_point() + ylim(c(0.4,1)) +
        ylab("Treatment Efficacy by AL") +
        xlab("Distance from Wild Type") +
        ggpubr::theme_pubclean() +
        theme(axis.line = element_line())

cowplot::plot_grid(g1,g2,g3, ncol = 3)

## alternative

df_grid <- data.frame(x = igraph::layout.reingold.tilford(g)[,1], y = igraph::layout.reingold.tilford(g)[,2],
                      dhappq = 1-drug_table$DHAPPQ,
                      al = 1-drug_table$AL,
                      asaq = 1-drug_table$ASAQ)

for (i in 0:6) {
        xs <- df_grid$x[df_grid$y == i]
        new_xs <- seq_along(xs)
        new_xs <- new_xs - median(new_xs)
        df_grid$x[df_grid$y == i] <- new_xs
}

dhappq <- ggplot(df_grid, aes(x,6-y,fill=1-dhappq)) + geom_tile(width = 1, color = "black") +
        scale_x_continuous("",expand = c(0,0)) +
        scale_y_reverse("Mutations from Wild Type",expand = c(0,0), position = "right") +
        scale_fill_gradientn("DHAPPQ 28-day \ntreatment \nsuccess",colours = terrain.colors(10)) +
        theme(panel.background = element_blank(),
              axis.ticks.x = element_blank(), axis.text.x = element_blank(),
              legend.position = "bottom") +
        coord_fixed()

al <- ggplot(df_grid, aes(x,6-y,fill=1-al)) + geom_tile(width = 1, color = "black") +
        scale_x_continuous("",expand = c(0,0)) +
        scale_y_reverse("Mutations from Wild Type",expand = c(0,0), position = "right") +
scale_fill_gradientn("AL 28-day \ntreatment \nsuccess",colours = terrain.colors(10)) +
        theme(panel.background = element_blank(),
              axis.ticks.x = element_blank(), axis.text.x = element_blank(),
              legend.position = "bottom") +
        coord_fixed()


#Lines
rayshader::plot_gg(dhappq, multicore = TRUE, raytrace = FALSE, width = 7, height = 4,
                   scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30)

rayshader::plot_gg(al, multicore = TRUE, raytrace = FALSE, width = 7, height = 4,
                   scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30)

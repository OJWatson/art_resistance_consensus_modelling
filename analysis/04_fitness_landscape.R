## Load in the drug by genotype table
load("/home/oj/GoogleDrive/AcademicWork/Imperial/git/magenta/data/drug_table.rda")

# Sort it into format to make a network
dt_df <- data.table::rbindlist(lapply(lapply(strsplit(drug_table$Genotype, ""),t),as.data.frame))
names(dt_df) <- c("pfcrt76T", "pfmdr186Y", "pfmdr1184F", "pfmdr1CNV", "pfk13580Y", "pfpm2CNV")
dt_df$Fitness <- as.numeric(drug_table$DHAPPQ)
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
colfunc <- colorRampPalette(colors = brewer.pal(11, "RdBu")[2:10])
brs_dhappq <- seq(min(1-drug_table$DHAPPQ), max(1-drug_table$DHAPPQ), length.out = 50)
brs_al <- seq(min(1-drug_table$AL), max(1-drug_table$AL), length.out = 50)
brs_asaq <- seq(min(1-drug_table$DHAPPQ), max(1-drug_table$ASAQ), length.out = 50)
cols <- rev(colfunc(length(brs_dhappq)))

# Create the DHA_PPQ graph
svg("analysis/plots/dha_ppq_network.svg", width = 18,height = 8)

plot(g, layout = igraph::layout.reingold.tilford(g),
     vertex.label = igraph::vertex_attr(g, "genotype"),
     vertex.color = cols[unlist(lapply(1-drug_table$DHAPPQ, function(x) {which.min(abs(x - brs_dhappq))}))],
     vertex.frame.color = c("white", "black")[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     vertex.frame.cex = c(1,2)[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     vertex.size = 6.5,
     edge.width = 1,
     edge.arrow.width = 0.1,
     vertex.size = 5,
     edge.arrow.size = 0.2,
     vertex.label.color = "black",
     vertex.size2 = 4,
     vertex.label.cex = 0.85,
     asp = 0.275,
     margin = -0.05)

image.plot(legend.only=T, zlim=range(brs_dhappq), col=cols,
           legend.lab = "\n\nProbability of 28-Day DHA-PPQ Treatment Failure",
           legend.mar = 3, horizontal = TRUE)
dev.off()

# Create the AL graph
svg("analysis/plots/al_network.svg", width = 18,height = 8)

plot(g, layout = igraph::layout.reingold.tilford(g),
     vertex.label = igraph::vertex_attr(g, "genotype"),
     vertex.color = cols[unlist(lapply(1-drug_table$AL, function(x) {which.min(abs(x - brs_al))}))],
     vertex.frame.color = c("white", "black")[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     vertex.frame.cex = c(1,2)[as.numeric(substr(igraph::vertex_attr(g, "genotype"), 6, 6)=="Y")+1],
     vertex.size = 6.5,
     edge.width = 1,
     edge.arrow.width = 0.1,
     vertex.size = 5,
     edge.arrow.size = 0.2,
     vertex.label.color = "black",
     vertex.size2 = 4,
     vertex.label.cex = 0.85,
     asp = 0.275,
     margin = -0.05)

image.plot(legend.only=T, zlim=range(brs_dhappq), col=cols,
           legend.lab = "\n\nProbability of 28-Day AL Treatment Failure",
           legend.mar = 3, horizontal = TRUE)
dev.off()

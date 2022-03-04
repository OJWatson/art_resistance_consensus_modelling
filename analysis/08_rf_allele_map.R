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
dt_df <- cbind(dt_df, !dt_df)
names(dt_df)[7:12] <- c("pfcrt76", "pfmdr186", "pfmdr1184", "pfmdr1", "pfk13580", "pfpm2")
dt_df <- dt_df %>% mutate_all(as.integer)

dt_df$tf_dhappq <- 1 - drug_table$DHAPPQ
dt_df$tf_al <- 1 - drug_table$AL
dt_df$tf_asaq <- 1 - drug_table$ASAQ


rf_dhappq <- randomForest::randomForest(
  tf_dhappq ~ .,
  data = dt_df[,c(1:12,13)], importance = TRUE
)

rf_al <- randomForest::randomForest(
  tf_al ~ .,
  data = dt_df[,c(1:12,14)], importance = TRUE
)

rf_asaq <- randomForest::randomForest(
  tf_asaq ~ .,
  data = dt_df[,c(1:12,15)], importance = TRUE
)

dt_df$dha_ppq_pred <- predict(rf_dhappq, newdata=dt_df[,1:12])
dt_df$al_pred <- predict(rf_al, newdata=dt_df[,1:12])
dt_df$asaq_pred <- predict(rf_asaq, newdata=dt_df[,1:12])
dt_df$muts <- apply(dt_df[,1:6], 1, sum)

g1 <- ggplot(dt_df, aes(muts, dha_ppq_pred, color = interaction(pfpm2CNV, pfk13580Y))) + geom_jitter(width=0.05) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limits = c(0,0.6)) + theme_bw() +
  scale_color_discrete(labels = c("Non-significant", "plasmepsin 2,3 CNV", "kelch13 580Y", "kelch13 580Y + plasmepsin 2,3 CNV"),
                       name = "")+
  xlab("Number of mutations from ancestral parasite") +
  ylab("Predicted Probability of treatment failure with DHA-PPQ") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA)) +
  guides(color=guide_legend(nrow=3,byrow=TRUE))


dt_df$AL_classifier <- "Non Significant"
dt_df$AL_classifier[dt_df$pfmdr1CNV == 0] <- "No pfmdr1 CNV"
dt_df$AL_classifier[dt_df$pfk13580Y == 1] <- "kelch13 580Y"
dt_df$AL_classifier[dt_df$pfk13580Y == 1 & dt_df$pfmdr1CNV == 0] <- "kelch13 580Y + No pfmdr1 CNV"
dt_df$AL_classifier[dt_df$pfk13580Y == 1 & dt_df$pfmdr186 == 1 & dt_df$pfmdr1184F == 1 & dt_df$pfmdr1CNV == 1] <- "kelch13 580Y + pfmdr1 86N 184F + pfmdr1 CNV"
dt_df$AL_classifier <- factor(
  dt_df$AL_classifier,
  levels = c("Non Significant", "No pfmdr1 CNV", "kelch13 580Y", "kelch13 580Y + No pfmdr1 CNV", "kelch13 580Y + pfmdr1 86N 184F + pfmdr1 CNV"))

g2 <- ggplot(dt_df, aes(muts, al_pred, color = AL_classifier)) + geom_jitter(width=0.05) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4), limits = c(0,0.45)) + theme_bw() +
  scale_color_discrete(name = "") +
  xlab("Number of mutations from ancestral parasite") +
  ylab("Predicted Probability of treatment failure with AL") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA)) +
  guides(color=guide_legend(nrow=3,byrow=TRUE))


g3 <- ggplot(dt_df, aes(muts, asaq_pred, color = interaction(pfmdr186Y, pfk13580Y))) + geom_jitter(width=0.05) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), limits = c(0,0.3)) + theme_bw() +
  scale_color_discrete(labels = c("Non-significant", "pfmdr1 86Y", "kelch13 580Y", "kelch13 580Y + pfmdr1 86Y"),
                       name = "")+
  xlab("Number of mutations from ancestral parasite") +
  ylab("Predicted Probability of treatment failure with ASAQ") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA)) +
  guides(color=guide_legend(nrow=3,byrow=TRUE))


varimps <- left_join(
  left_join(
  data.frame(allele = rownames(rf_dhappq$importance),
                      imp_dhappq = rf_dhappq$importance[,2]),
  data.frame(allele = rownames(rf_al$importance),
             imp_al = rf_al$importance[,2]), by = "allele"),
  data.frame(allele = rownames(rf_asaq$importance),
             imp_asaq = rf_asaq$importance[,2]), by = "allele")

varimps$allele <- c("pfcrt 76T", "pfmdr1 86Y", "pfmdr1 184F",
                    "pfmdr1 CNV", "pfk13 580Y", "pfpm2,3 CNV",
                    "pfcrt 76K", "pfmdr1 86N", "pfmdr1 184Y",
                    "pfmdr1 single copy", "pfk13 580C", "pfpm2,3 single copy")

v1 <- ggplot(varimps, aes(allele, imp_dhappq)) + geom_bar(stat = "identity") +
  xlab("Allele") +
  ylab("Random Forest Variable Importance for DHAPPQ") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA),axis.text.x = element_text(angle = 45, hjust = 1))

v2 <- ggplot(varimps, aes(allele, imp_al)) + geom_bar(stat = "identity") +
  xlab("Allele") +
  ylab("Random Forest Variable Importance for AL") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA),axis.text.x = element_text(angle = 45, hjust = 1))

v3 <- ggplot(varimps, aes(allele, imp_asaq)) + geom_bar(stat = "identity") +
  xlab("Allele") +
  ylab("Random Forest Variable Importance for ASAQ") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.key = element_rect(fill = NA),axis.text.x = element_text(angle = 45, hjust = 1))

## Random forest performance

df_res <- data.frame("prediction" = c(rf_dhappq$predicted, rf_asaq$predicted, rf_al$predicted),
           "observed" = c(dt_df$tf_dhappq, dt_df$tf_asaq, dt_df$tf_al),
           "drug" = c(rep("DHA-PPQ", 64),rep("ASAQ", 64),rep("AL", 64)))


rf_gg <- ggplot(df_res, aes(observed, prediction, color = drug)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("Observed 28-Day Treatment Failure") +
  ylab("Predicted 28-Day Treatment Failure") +
  scale_color_viridis_d(name = "Drug") +
  theme_bw()

rf_gg_2 <- df_res %>% group_by(drug) %>%
  summarise(`1 - Correlation` = 1-round(cor(prediction, observed), 4),
            MAE = round(mean(abs(prediction - observed)),4),
            RMSE = round(sqrt(mean((prediction - observed)^2)),4)) %>%
  pivot_longer(`1 - Correlation`:RMSE) %>%
  ggplot(aes(name, value, fill = drug)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  xlab("Model Metric") +
  ylab("Model Performance") +
  scale_fill_viridis_d(name = "Drug") +
  theme_bw()


comb <- cowplot::plot_grid(g1,g2,g3,NA, NA, NA, v1,v2,v3, ncol = 3, labels = c("a","b","c","","","","d","e","f"), rel_heights = c(1,0.1,1))
comb2 <- cowplot::plot_grid(rf_gg, rf_gg_2, labels = c("h", "i"), ncol = 2)

save_figs("rf_allele", cowplot::plot_grid(comb, NA, comb2, ncol = 1, rel_heights = c(0.75, 0.025, 0.25)), height = 16, width = 16)

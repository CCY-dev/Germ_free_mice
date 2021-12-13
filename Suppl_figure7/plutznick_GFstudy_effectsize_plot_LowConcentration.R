rm(list = ls())
setwd("~/Documents/Lab/Ellen_gremFree/analysis/")
library(tidyverse)
library(ggpubr)
library(orddom)

GF_metabolite <- read.table("../rawdata_fromDropbox/Low_conc_metabolite_comparison_data/full_wGFandCOL_imputedlowcondata_matching_15102020_edited.csv", header = T,
                            sep = ",", row.names = NULL, stringsAsFactors = F)
GF_metabolite_reduced <- GF_metabolite %>%
  mutate(Study = "Berlin") %>%
  dplyr::select(ID, Study, treatment, microbiome, everything())
GF_metabolite_reduced_long <- GF_metabolite_reduced %>%
  pivot_longer(5:ncol(GF_metabolite_reduced))


Plutz_metabolite <- read.table("../rawdata_fromDropbox/Low_conc_metabolite_comparison_data/Plutzdata_simplified_with_matching_15102020_edited.csv", header = T,
                               sep = ",", row.names = NULL, stringsAsFactors = F)
Plutz_reduced <- Plutz_metabolite %>%
  mutate(Study = "Plutznick") %>%
  dplyr::select(ID, Study, treatment, microbiome, everything())
Plutz_reduced_long <- Plutz_reduced %>%
  pivot_longer(5:ncol(Plutz_reduced))

metabolites <- unique(GF_metabolite_reduced_long$name)


######Prepare Ellen's data, cliff's delta (Ang to Sham)
GF_effectsize_GF <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(GF_effectsize_GF) <- c("Metabolite", "GF_Ang2_to_sham_effect_size", "Study")
GF_effectsize_GF$Metabolite <- metabolites
GF_effectsize_GF$Study <- "Berlin"
for (i in 1:length(metabolites)) {
  GF_sub_GF <- GF_metabolite_reduced_long %>%
    filter(name == metabolites[i], microbiome == "GF")
  GF_sub_GF_Ang <- GF_sub_GF %>% filter(treatment == "AngII")
  GF_sub_GF_Sham <- GF_sub_GF %>% filter(treatment == "Sham")
  GF_effectsize_GF[i, 2] <- dmes(GF_sub_GF_Ang$value, GF_sub_GF_Sham$value)$dc
}

GF_effectsize_COL <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(GF_effectsize_COL) <- c("Metabolite", "COL_Ang2_to_sham_effect_size", "Study")
GF_effectsize_COL$Metabolite <- metabolites
GF_effectsize_COL$Study <- "Berlin"
for (i in 1:length(metabolites)) {
  GF_sub_COL <- GF_metabolite_reduced_long %>%
    filter(name == metabolites[i], microbiome == "COL")
  GF_sub_COL_Ang <- GF_sub_COL %>% filter(treatment == "AngII")
  GF_sub_COL_Sham <- GF_sub_COL %>% filter(treatment == "Sham")
  GF_effectsize_COL[i, 2] <- dmes(GF_sub_COL_Ang$value, GF_sub_COL_Sham$value)$dc
}

#Prepare Plutznick's male data
PL_effectsize_GF <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(PL_effectsize_GF) <- c("Metabolite", "GF_Ang2_to_sham_effect_size", "Study")
PL_effectsize_GF$Metabolite <- metabolites
PL_effectsize_GF$Study <- "Plutznick"
for (i in 1:length(metabolites)) {
  PL_sub_GF <- Plutz_reduced_long %>%
    filter(name == metabolites[i], microbiome == "GF")
  PL_sub_GF_Ang <- PL_sub_GF %>% filter(treatment == "AngII")
  PL_sub_GF_Sham <- PL_sub_GF %>% filter(treatment == "Sham")
  PL_effectsize_GF[i, 2] <- dmes(PL_sub_GF_Ang$value, PL_sub_GF_Sham$value)$dc
}


PL_effectsize_COL <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(PL_effectsize_COL) <- c("Metabolite", "COL_An2_to_sham_effect_size", "Study")
PL_effectsize_COL$Metabolite <- metabolites
PL_effectsize_COL$Microbiome <- "COL"
PL_effectsize_COL$Study <- "Plutznick"
for (i in 1:length(metabolites)) {
  PL_sub_COL <- Plutz_reduced_long %>%
    filter(name == metabolites[i], microbiome == "COL")
  PL_sub_COL_Ang <- PL_sub_COL %>% filter(treatment == "AngII")
  PL_sub_COL_Sham <- PL_sub_COL %>% filter(treatment == "Sham")
  PL_effectsize_COL[i, 2] <- dmes(PL_sub_COL_Ang$value, PL_sub_COL_Sham$value)$dc
}

# Combine the tables together
Ellen_effectsize <- as.data.frame(cbind(GF_effectsize_GF, GF_effectsize_COL$COL_Ang2_to_sham_effect_size)) %>%
  dplyr::select(1, 2, 4, 3)
colnames(Ellen_effectsize)[3] <- "COL_Ang2_to_Sham_effect_size"
#write.table(Ellen_effectsize, "Ellen_effectsize.txt", sep = "\t")

PL_effectsize <- as.data.frame(cbind(PL_effectsize_GF, PL_effectsize_COL$COL_An2_to_sham_effect_size)) %>%
  dplyr::select(1, 2, 4, 3)
colnames(PL_effectsize)[3] <- "COL_Ang2_to_Sham_effect_size"
#write.table(PL_effectsize, "PL_effectsize.txt", sep = "\t")

All_effectsize <- rbind(Ellen_effectsize, PL_effectsize)

# Plotting

#pdf("GF_COL_effectsize_pairedLine_lowConc.pdf", width = 11, height = 6)
ggplot(data = All_effectsize, aes(x = GF_Ang2_to_sham_effect_size, y = COL_Ang2_to_Sham_effect_size))+
  geom_point(aes(color = Study), position = position_nudge(x = 0),
             alpha = 0.8, size = 3) +
  theme_classic() + ggtitle("GF v.s. COL effect size") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.y = element_text(face = "italic")) +
  xlab(label = "Ang2 to sham effect size in germ free mice ") +
  ylab(label = "Ang2 to sham effect size in colonized mice ") +
  coord_fixed(ratio = 1) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  geom_density_2d(aes(color = Study), size = 0.35)
#dev.off()

# Calculate the distance between metabolites of 2 studies on x and y axis
distance <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(distance) <- c("Metabolite", "Dist_germ_free", "Dist_colonized")
distance$Metabolite <- metabolites
for (i in 1:length(metabolites)) {
  sub <- All_effectsize %>%
    filter(Metabolite == metabolites[i])
  distance[i, 2] <- abs(sub[2, 2] - sub[1, 2])
  distance[i, 3] <- abs(sub[2, 3] - sub[1, 3])
}
distance_long <- distance %>%
  pivot_longer(2:3)
wilcox.test(distance$Dist_germ_free, distance$Dist_colonized, paired = T)

library(ggpubr)
#pdf("Difference_in_metabolite_effect_sizes_Berlin_Plustznick_LowConc.pdf",
#    width = 6, height = 7)
ggplot(distance_long, aes(name, value)) +
  geom_violin(aes(fill = name)) +
  geom_boxplot(width = 0.2) + # Combine with box plot to add median and quartiles
  ggtitle("Difference in effect sizes") +
  theme_light() + xlab("Group") + ylab("Absolute distance") +
  stat_compare_means(inherit.aes = TRUE, method = "wilcox.test", paired = T,
                     aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  expand_limits(y = c(0, 2.3)) +
  theme(legend.position = "none") #+ # Remove legend
  #geom_line(aes(name, value, group = Metabolite))
#dev.off()













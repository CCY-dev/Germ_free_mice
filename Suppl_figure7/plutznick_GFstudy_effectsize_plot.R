rm(list = ls())
setwd("~/Documents/Lab/Ellen_gremFree/analysis/")
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(orddom)

GF_metabolite <- read.table("../rawdata_fromDropbox/metabolite comparison data/Our data/peruMmetabolitedata_forcomp_18062020.csv", header = T,
                          sep = ",", row.names = NULL, stringsAsFactors = F)

# Drop the ones whose MATCH_NAME = NA
# Also replace na in value column with 0
GF_metabolite_reduced <- GF_metabolite %>%
  drop_na(MATCH_NAME) %>%
  arrange(MATCH_NAME) %>%
  dplyr::select(-shortform_biocrates) %>%
  mutate(Study = "Ellen") %>%
  dplyr::select(id, Study, treatment, microbiome, MATCH_NAME, value) %>%
  mutate_all(~replace(., is.na(.), 0))

Plutz_metabolite <- read_xlsx(path ="../rawdata_fromDropbox/metabolite comparison data/Plutznick/VolnormImput_reduced.xlsx",
                              col_names = T)
Plutz_long <- Plutz_metabolite %>%
  gather(key = "MATCH_NAME", value = "value", c(5:57)) %>%
  arrange(MATCH_NAME) %>%
  mutate(Study = "Plutznick") %>%
  dplyr::select(id, Study, treatment, microbiome, MATCH_NAME, value, Gender)

# Select only the males in Plutz study!!!!!
Plutz_long_M <- Plutz_long %>%
  filter(Gender == "M") %>%
  dplyr::select(-Gender)

metabolites <- unique(Plutz_long_M$MATCH_NAME)

#Prepare Ellen's data, cliff's delta (Ang to Sham)
GF_effectsize_GF <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(GF_effectsize_GF) <- c("Metabolite", "GF_Ang2_to_sham_effect_size", "Study")
GF_effectsize_GF$Metabolite <- metabolites
GF_effectsize_GF$Study <- "Ellen"
for (i in 1:length(metabolites)) {
  GF_sub_GF <- GF_metabolite_reduced %>%
    filter(MATCH_NAME == metabolites[i], microbiome == "GF")
  GF_sub_GF_Ang <- GF_sub_GF %>% filter(treatment == "AngII")
  GF_sub_GF_Sham <- GF_sub_GF %>% filter(treatment == "Sham")
  GF_effectsize_GF[i, 2] <- dmes(GF_sub_GF_Ang$value, GF_sub_GF_Sham$value)$dc
}

GF_effectsize_COL <- as.data.frame(matrix(NA, nrow = length(metabolites), ncol = 3))
colnames(GF_effectsize_COL) <- c("Metabolite", "COL_Ang2_to_sham_effect_size", "Study")
GF_effectsize_COL$Metabolite <- metabolites
GF_effectsize_COL$Study <- "Ellen"
for (i in 1:length(metabolites)) {
  GF_sub_COL <- GF_metabolite_reduced %>%
    filter(MATCH_NAME == metabolites[i], microbiome == "COL")
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
  PL_sub_GF <- Plutz_long_M %>%
    filter(MATCH_NAME == metabolites[i], microbiome == "GF")
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
  PL_sub_COL <- Plutz_long_M %>%
    filter(MATCH_NAME == metabolites[i], microbiome == "COL")
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

#pdf("GF_COL_effectsize_pairedLine_male.pdf", width = 11, height = 6)
ggplot(data = All_effectsize, aes(x = GF_Ang2_to_sham_effect_size, y = COL_Ang2_to_Sham_effect_size))+
  geom_point(aes(color = Study), position = position_nudge(x = 0),
             alpha = 0.8, size = 3) +
  theme_classic() + ggtitle("GF v.s. COL effect size (male)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.y = element_text(face = "italic")) +
  xlab(label = "Ang2 to sham effect size in germ free mice ") +
  ylab(label = "Ang2 to sham effect size in colonized mice ") +
  coord_fixed(ratio = 1) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  #geom_line(aes(x = GF_Ang2_to_sham_effect_size,
  #              y = COL_Ang2_to_Sham_effect_size, group = Metabolite),
  #          color = "lightgrey") +
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
#pdf("Difference_in_metabolite_effect_sizes_Ellen_Plustznick_male.pdf",
#    width = 6, height = 7)
ggplot(distance_long, aes(name, value)) +
  geom_violin(aes(fill = name)) +
  geom_boxplot(width = 0.2) + # Combine with box plot to add median and quartiles
  ggtitle("Difference in effect sizes") +
  theme_light() + xlab("Group") + ylab("Absolute distance") +
  stat_compare_means(inherit.aes = TRUE, method = "wilcox.test", paired = T,
                     aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  expand_limits(y = c(0, 2.3)) +
  theme(legend.position = "none")# + # Remove legend
  #geom_line(aes(name, value, group = Metabolite))
#dev.off()

########### Do PCoA on the effect size to see if GF cluster more than COL
#library(ade4)
#library(vegan)
#library(ape)
#dist <- vegdist(distance_long[ , 3],  method = "euclidean", na.rm = T)
#all.pcoa <- cmdscale(dist, k = (nrow(distance_long[ , c(3)])-1),
#                     eig=TRUE, add = F)
#in_data <- as.data.frame(vegan::scores(all.pcoa, choices = c(1,2)))
#in_data$Group <- distance_long$name


#pdf("PCoA_effectsize_male.pdf", width = 8, height= 6 )
#ggplot(data = in_data, aes(x = Dim1, y = Dim2, color = Group)) +
#  geom_point(size = 3, alpha = 0.8) + theme_classic() + ggtitle("PCoA effect size male (euclidean)") +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  stat_ellipse(aes(x=Dim1,y=Dim2,color = Group),level = 0.50)

#dev.off()


############# Adonis test
#distance_long_new <- distance_long
#colnames(distance_long_new)[2:3] <- c("Group", "distance")
#adonis(data = distance_long_new, formula = distance ~ Group)

############## Plot NMDS
#all.mds <- metaMDS(distance_long_new[ , 3])
#data.scores <- as.data.frame(scores(all.mds))
#data.scores$site <- rownames(data.scores)
#data.scores$Group <- distance_long$name


#pdf("NMDS_TS42.pdf", width = 8, height= 6 )
#ggplot(data = data.scores) +
#  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=Group),level = 0.95) +
#  geom_point(aes(x=NMDS1,y=NMDS2,colour=Group),size=1.5, alpha = 0.8) +
#  theme_classic() + ggtitle("NMDS TS4") #+
#geom_text(aes(x = NMDS1, y = NMDS2, label = Sample), nudge_x = 0.07) +
#geom_line(aes(x = NMDS1, y = NMDS2, group = Sample), color = "darkgrey")
#dev.off()



# Plotting
#pdf("GF_COL_effectsize_label.pdf", width = 11, height = 6)
#ggplot(data = All_effectsize, aes(x = GF_Ang2_to_sham_effect_size, y = COL_Ang2_to_Sham_effect_size))+
#  geom_point(aes(color = Study), position = position_nudge(x = 0),
#             alpha = 0.7, size = 2) +
#  theme_classic() + ggtitle("GF v.s. COL effect size") +
#  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
#  theme(axis.text.y = element_text(face = "italic")) +
#  xlab(label = "Ang2 to sham effect size in germ free mice ") +
#  ylab(label = "Ang2 to sham effect size in colonized mice ") +
#  coord_fixed(ratio = 1) +
#  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
#  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
#  geom_text_repel(aes(label = Metabolite), size = 2.5)
#dev.off()











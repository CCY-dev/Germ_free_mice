rm(list = ls())
setwd("~/Documents/Lab/Ellen_gremFree/analysis/")
library(tidyverse)
library(phytools)
library(vegan)
library(ade4)
library(readxl)

# Read in Ellen_GF study meta data
metadata <- read.table("/Users/Jessica/Documents/Lab/Ellen_gremFree/rawdata_fromDropbox/metadata.txt", header = T,
                       sep = "\t", row.names = NULL, stringsAsFactors = F)

# Read in Ellen_GF study microbes genus and meta data
GFmicrobe <- read.table("/Users/Jessica/Documents/Lab/Ellen_gremFree/rawdata_fromDropbox/GF_study_shotgun_result/output.genus.rrarefied_to_1130.000000_n_0.tsv", sep = "\t")
GFmicrobe <- as.data.frame(t(GFmicrobe))
names(GFmicrobe) <- GFmicrobe %>% slice(1) %>% unlist()
GFmicrobe <- GFmicrobe %>% slice(-1)
colnames(GFmicrobe)[1] <- "Sample"
GFmicrobe_all <- inner_join(metadata, GFmicrobe, by = "Sample")
GFmicrobe_feces <- GFmicrobe_all %>% filter(Origin == "feces")
GFmicrobe_feces <- GFmicrobe_feces %>% dplyr::select(-c(2, 3, 4, 6))
GFmicrobe_feces <- GFmicrobe_feces %>%
  mutate(Sex = "M") %>%
  dplyr::select(Sample, Treatment, Sex, everything())
# Remove unknown/unassigned
GFmicrobe_feces$Unassigned <- NULL
GFmicrobe_feces$Study <- "GF"

# Read in Plutznick microbes genus
Plutznick <- read.table("/Users/Jessica/Documents/Lab/Ellen_gremFree/rawdata_fromDropbox/Plutznick_16S/Genus_rtk.txtrarefied_to_169838.000000_n_0_onlyGenusName.tsv", sep = "\t")
Plutznick <- as.data.frame(t(Plutznick))
names(Plutznick) <- Plutznick %>% slice(1) %>% unlist()
Plutznick <- Plutznick %>% slice(-1)
colnames(Plutznick)[1] <- "Sample"
Names <- as.data.frame(str_split_fixed(Plutznick$Sample, pattern = "_", n = 3))
colnames(Names) <- c("Treatment", "Sample", "Sex")
Plutznick_all <- as.data.frame(cbind(Names, Plutznick[ , -1])) %>%
  dplyr::select(Sample, Treatment, Sex, everything())
Plutznick_all$Study <- "Plutznick"

# Bind Ellen_GF and Plutznick together, and unify their data
GFmicrobe_feces_t <- as.data.frame(t(GFmicrobe_feces)) %>% 
  rownames_to_column("X")
Plutznick_all_t <- as.data.frame(t(Plutznick_all)) %>% 
  rownames_to_column("X")
Merge_t <- as.data.frame(inner_join(GFmicrobe_feces_t, Plutznick_all_t, by = "X"))

# Write table and do rarefaction on VM
#write.table(Merge_t, "GF_Plutznick_genus_merge.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# Read in the rarefied table
rarefied <- read.table("GF_Plutznick_genus_merge_rtkrarefied_to_181.000000_n_0.tsv", header = F,
                       sep = "\t", row.names = NULL, stringsAsFactors = F)
rarefied <- as.data.frame(apply(rarefied, 2, as.character))
colnames(rarefied) <- colnames(Merge_t)
rarefied_allinfo <- dplyr::bind_rows(rarefied, Merge_t[c(2, 3, 23), ])

Merge <- as.data.frame(t(rarefied_allinfo))
names(Merge) <- Merge %>% slice(1) %>% unlist()
Merge <- Merge %>% 
  slice(-1) %>% 
  dplyr::select(Sample, Study, Treatment, Sex, everything())

Merge_value <- Merge[ , -c(1:4)]
Merge_value <- apply(Merge_value, 2, as.character)
Merge_value <- apply(Merge_value, 2, as.numeric)

Treatment_original <- Merge$Treatment
Treatment_1 <- ifelse(Treatment_original == "Ang" | Treatment_original == "Ang2",
                      yes = "AngII", no = "Sham")
Merge$Treatment_correct <- Treatment_1
Merge <- Merge %>%
  dplyr::select(Sample, Study, Treatment, Treatment_correct, Sex, everything())


############ Metabolome PCoA on Cage effect for ALL mice! ############
dist <- vegdist(Merge_value,  method = "bray", na.rm = T)
all.pcoa <- cmdscale(dist, k = (nrow(Merge_value)-1), eig=TRUE)
#get the % of variation explained by each of the axis' using the eigen value 
mds.var.per <- round(all.pcoa$eig/sum(all.pcoa$eig)*100, 1)
in_data <- as.data.frame(vegan::scores(all.pcoa, choices = c(1,2)))
in_data$Treatment <- Merge$Treatment_correct
in_data$Study <- Merge$Study
in_data$Sex<- Merge$Sex

#pdf("PCoA_Plutznick_GF_mcirobiome_bray_genus_rarefied3.pdf", width = 8, height= 6 )
ggplot(data = in_data, aes(x = Dim1, y = Dim2, color = Study, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.7) + theme_classic() +
  ggtitle("PCoA of microbiome at genus level (bray)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(aes(color = Study, lty =  Treatment)) +
  xlab(paste("MDS1-", mds.var.per[1], "%", sep = "")) + 
  ylab(paste("MDS2-", mds.var.per[2], "%", sep = "")) 
dev.off()





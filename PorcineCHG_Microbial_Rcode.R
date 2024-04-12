# R code for microbial data:
# PMAxx optomization experiment 
# ex vivo skin CHG experiment, CFU and viability-qPCR bioburden, and 16S analysis 

library(phyloseq)
library(broom)
library(ggrepel)
library(ggplot2)
library(microbiome)
library(microbiomeutilities) # for the core microbiome analysis 
library(metagMisc)
library(tidyverse)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(decontam)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(Maaslin2)
library(reshape2)
library(pheatmap)
library(speedyseq)

setwd("/Users/liztown/Documents/KalanLab/Papers/PigSkinCHG/Rcode_FINAL")

# Supplemental Figure 2 - PMAxx optomization experiment
PCR <- read.csv("./PMAxxOptAllResults.csv")
PCR <- subset.data.frame(PCR, LHKC != "Heated") # removing the 90C for 5 minutes boil condition 

ggplot(data=PCR, aes(x=reorder(LHKC, LHKCorder)), color = factor(LHKC)) + 
  geom_boxplot(aes(y = MeanPercentLive,color = factor(LHKC), fill = factor(LHKC)), size = 1.3, alpha = 0.85)+
  geom_point(aes(y =MeanPercentLive, shape = Pig, color = factor(LHKC)),  size = 2.5)+
  facet_wrap(~PMAxxConcentration, scales = "free_x", ncol = 4)+
  scale_y_continuous(limits = c(0,60), breaks = c(0, 10, 20, 30, 40, 50, 60)) + 
  theme_light() +
  scale_fill_manual(values = c("#F47439","#FDB632","#023B58" ))+ # TYOB
  scale_color_manual(values = c("#C1450B","#CA8402", "#011B28"))+
  ggtitle("Viability PCR Pig Skin Optomization Part 1; % Live") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")


# Figure `1: Ex vivo porcine CHG experiment microbial bioburden `
# Figure 1C - Colony forming units 
CFU <- read.csv("./Aim2CFU_ALL.csv")
CFU$PigTxTime <- paste(CFU$Pig, CFU$TxTime)
CFU<- CFU %>% group_by(PigTxTime) %>%  mutate(MeanLogCFU= mean(LogCFU))

CFU$Treatment<- factor(CFU$Treatment, levels = c("Water", "Chloraprep", "Full"))

ggplot(data=CFU, aes(x=Timepoint, y=(MeanLogCFU), color = Treatment)) + 
  stat_summary(fun = "median", geom = "line", linetype = 2, linewidth = 1)  +
  geom_boxplot(aes(x=Timepoint, y= MeanLogCFU, group = TxTime, fill = Treatment, color = Treatment), alpha = 0.85, size = 1, width = 6, position = position_dodge2(1))+
  scale_y_continuous( breaks = c(0, 2, 4, 6, 8)) + 
  scale_fill_manual(values = c("#027979","#FDB632","#F47439"))+
  scale_color_manual(values = c("#012828","#B67702","#9A3709"))+
  scale_x_continuous(breaks = c(-12,0, 6, 12, 24, 48)) + 
  geom_hline(yintercept=c(0.88), size = 1, linetype="dotted") +
  theme_light() +
  ggtitle("Aim 2; CFU results. all 3 pigs") +
  ylab("Log10(CFU)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))

#Figure 1 D-E- viability PCR - viable and total bioburden ex vivo CHG experiment 
PCR <- read.csv("./Aim2-qPCRallPigs.csv")
PCR <- as.data.frame(PCR)
# Negative controls  == 174, 105,and 110 bacteria, for Pig 1,2,and3 qPCR respectively. average == 130
PCR <- subset.data.frame(PCR,SampleName != "NegEX" )
PCR$Treatment2 <-factor(PCR$Treatment2, levels = c("Baseline", "Water", "Chloroprep", "CHG Bath (x2) + Chloroprep"))

ggplot(data=PCR, aes(x=Time, color = Treatment2)) + 
  geom_boxplot(aes(x = Time, y = log10(ViableMeanBacteria), group = Group2, fill = Treatment2), size = 1, width = 6, position = position_dodge(-1.7))+
  stat_summary(aes(y = log10(ViableMeanBacteria)), fun = "mean", geom = "line", size = 1., linetype = 2, position = position_dodge(-1.7)) +
  scale_linetype_manual("", values=c("dashed"))+
  scale_fill_manual(values = c("#CCCCCC","#027979","#FDB632","#F47439"))+
  scale_color_manual(values = c("#777777","#012828","#B67702","#9A3709"))+
  scale_x_continuous(breaks = c(-12, 0, 6, 12, 24, 48), minor_breaks = NULL) + 
  geom_hline(yintercept=c(log10(130)), size = 1, linetype="dotted") +
  theme_light() +
  ggtitle("Aim2 Experiment PCR Calculated Bacteria") +
  ylab("Calculated Bacteria Log10(Bacteria) - all 3 pigs combined") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data=PCR, aes(x=Time, color = Treatment2)) + 
  geom_boxplot(aes(x = Time, y = log10(TotalMeanBacteria), group = Group2, fill = Treatment2), size = 1, width = 6, position = position_dodge(1.7))+
  stat_summary(aes(y = log10(TotalMeanBacteria)), fun = "mean", geom = "line", size = 1.,position = position_dodge(1.7)) +
  scale_linetype_manual("", values=c("solid"))+
  scale_fill_manual(values = c("#CCCCCC","#027979","#FDB632","#F47439"))+
  scale_color_manual(values = c("#777777","#012828","#B67702","#9A3709"))+
  scale_x_continuous(breaks = c(-12, 0, 6, 12, 24, 48), minor_breaks = NULL) + 
  geom_hline(yintercept=c(log10(130)), size = 1, linetype="dotted") +
  theme_light() +
  ggtitle("Aim2 Experiment PCR Calculated Bacteria") +
  ylab("Calculated Bacteria Log10(Bacteria) - all 3 pigs combined") +
  theme(plot.title = element_text(hjust = 0.5))



# 16S analysis for all the pigs 
# ps.NAME = a phyloseq object --> this is importing the datasets and associated metadata 
ps.Aim2_P1 <- qza_to_phyloseq(features = "./LK16S005-table-Aim2-dada2.qza",
                           taxonomy = "./LK16S005-taxonomy-Aim2.qza",
                           metadata = "./LK16S005-Aim2_P1_SequencingManufest.txt")
ps.Aim2_P2 <-  qza_to_phyloseq(features = "./LK16S008-1to96-table-16S-dada2.qza",
                                taxonomy = "./LK16S008-1to96-taxonomy-Aim2P2.qza",
                                metadata = "./LK16S008_16S_1to96_manufest.txt")
ps.Aim2_P3 <-  qza_to_phyloseq(features = "./LK16S008-97plus-table-16S-dada2.qza",
                               taxonomy = "./LK16S008-97plus-taxonomy-Aim2P3.qza",
                               metadata = "./LK16S008_16S_97plus_manufest.txt")
ps.Aim2_P3re <- qza_to_phyloseq(features = "./LK16S011-table-16S-dada2.qza",
                                taxonomy = "./LK16S011-taxonomy-Aim2P3re.qza",
                                metadata = "./LK16S011_P3_16Sreruns_manufest.txt")
                               
# function to remove specific bad / contaminant OTUs/ASVs
remove_BadOTU = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

# Removing contaminants 
#ps.Decontam <- prune_taxa(!contamdf.prev$contaminant, ps.PSE) # REMOVING THE CONTAMINANTS identified via DECONTAM 
ps.Decontam <- subset_taxa(ps.Aim2_P1, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 

ps.Decontam <- subset_taxa(ps.Decontam, Genus !="0319-6G20") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Acidibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Anaeromyxobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Bradyrhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curtobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Dyadobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Enhydrobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Entotheonellaceae") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Exiguobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ferruginibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gaiella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gracilibacteria") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Herbaspirillum") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="LWQ8") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Marvinbryantia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="MB-A2-108") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylotenera") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Motilibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Nakamurella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="NS9_marine_group") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rathayibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rhodoferax") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sandaracinobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sediminibacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingopyxis") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="SWB02") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Tepidiphilus") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Thermomonas") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Undibacterium") 

ps.Decontam <- subset_samples(ps.Decontam, !sample_names(ps.Decontam) %in%  c("LK16S005-097", "LK16S005-098")) # removing the negative control
ps.P1 <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.P1# THE CLEANED FILE 

bad.OTU <- c("8e38d45e13733c913016c40b14a3bb87")  # prominant OTU in the negative controls 
Aim2_P3re_Clean <- remove_BadOTU(ps.Aim2_P3re, bad.OTU)
ps.Decontam <- subset_taxa(Aim2_P3re_Clean, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus!="Clostridium") # pmaxx negative extraction contaminant 

ps.Decontam <- subset_taxa(ps.Decontam, Genus !="0319-6G20") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Acidibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Anaeromyxobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Bradyrhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curtobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Dyadobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Enhydrobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Entotheonellaceae") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Exiguobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ferruginibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gaiella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gracilibacteria") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Herbaspirillum") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="LWQ8") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Marvinbryantia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="MB-A2-108") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylotenera") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Motilibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Nakamurella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="NS9_marine_group") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rathayibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rhodoferax") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sandaracinobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sediminibacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingopyxis") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="SWB02") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Tepidiphilus") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Thermomonas") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Undibacterium") 

ps.Decontam <- subset_samples(ps.Decontam, !sample_names(ps.Decontam) %in% c("LK16S008-199", "LK16S008-200")) # removing the negative control
ps.P3re <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.P3re@sam_data$Pig[ps.P3re@sam_data$Pig == "Pig3"] <- "P3"

bad.OTU <- c("8e38d45e13733c913016c40b14a3bb87")  # prominant OTU in the negative controls 
Aim2_P3_Clean <- remove_BadOTU(ps.Aim2_P3, bad.OTU)
ps.Decontam <- subset_taxa(Aim2_P3_Clean, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus!="Clostridium") # pmaxx negative extraction contaminant 

ps.Decontam <- subset_taxa(ps.Decontam, Genus !="0319-6G20") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Acidibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Anaeromyxobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Bradyrhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curtobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Dyadobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Enhydrobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Entotheonellaceae") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Exiguobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ferruginibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gaiella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gracilibacteria") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Herbaspirillum") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="LWQ8") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Marvinbryantia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="MB-A2-108") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylotenera") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Motilibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Nakamurella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="NS9_marine_group") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rathayibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rhodoferax") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sandaracinobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sediminibacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingopyxis") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="SWB02") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Tepidiphilus") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Thermomonas") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Undibacterium") 

ps.Decontam <- subset_samples(ps.Decontam, !sample_names(ps.Decontam) %in% c("LK16S008-199", "LK16S008-200")) # removing the negative control
ps.DeP3 <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.P3.1 <- subset_samples(ps.DeP3, Pig == "P3")
ps.P2.2<-subset_samples(ps.DeP3, Pig == "P2")
ps.P3<- merge_phyloseq(ps.P3.1, ps.P3re)
ps.P3 # FINAL phyloseq object

bad.OTU <- c("8e38d45e13733c913016c40b14a3bb87")  # prominant OTU in the negative controls 
Aim2_P2_Clean <- remove_BadOTU(ps.Aim2_P2, bad.OTU)
ps.Decontam <- subset_taxa(Aim2_P2_Clean, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus!="Clostridium") # pmaxx negative extraction contaminant 

ps.Decontam <- subset_taxa(ps.Decontam, Genus !="0319-6G20") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Acidibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Anaeromyxobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Bradyrhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curtobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Dyadobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Enhydrobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Entotheonellaceae") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Exiguobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ferruginibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gaiella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gracilibacteria") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Herbaspirillum") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="LWQ8") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Marvinbryantia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="MB-A2-108") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylotenera") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Motilibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Nakamurella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="NS9_marine_group") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rathayibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rhodoferax") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sandaracinobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sediminibacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingopyxis") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="SWB02") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Tepidiphilus") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Thermomonas") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Undibacterium") 

ps.DeP2 <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.P2 <- merge_phyloseq(ps.DeP2, ps.P2.2)

ps.Aim2 <- merge_phyloseq(ps.P1, ps.P2, ps.P3)
ps.Aim2@sam_data$Timepoint[ps.Aim2@sam_data$Timepoint == "Post Tx"] <- "0 hr"
ps.Aim2 <- subset_samples(ps.Aim2, Pig != "20 strain")
ps.Aim2@sam_data$Base.WaterCHG <- ps.Aim2@sam_data$WaterCHG
ps.Aim2@sam_data$Base.WaterCHG[ps.Aim2@sam_data$Timepoint == "Baseline"] <- "Baseline"
ps.Aim2@sam_data$Base.WaterCHG 

OutCV <- psmelt(ps.Aim2)
#write.csv(OutCV, "./Aim2_AllPigs_AbundanceTotalsOut.csv")

ps.Aim2_REl <- transform(ps.Aim2, "compositional")
Aim2.Relout <- psmelt(ps.Aim2_REl)
#write.csv(Aim2.Relout, "./Aim2_AllPigs_RELATIVEAbundanceTotalsOut.csv")

Aim2.Rel.df<-as.data.frame(Aim2.Relout)
Aim2.meta <- data.frame(sample_data(ps.Aim2))
Aim2.meta$PMAxxSample <- row.names(Aim2.meta)


# ALPHA ABUNDANCE (Not shown in paper)
ps.Aim2@sam_data$Timepoint <- factor(ps.Aim2@sam_data$Timepoint, levels = c("Baseline", "0 hr","6 hr","12 hr", "24 hr", "48 hr"))
ps.Aim2@sam_data$WaterCHG <- factor(ps.Aim2@sam_data$WaterCHG, levels = c("Water", "CHG", "Full"))
ps.Aim2@sam_data$PmaxxGroup <- factor(ps.Aim2@sam_data$PmaxxGroup, levels = c("Water","Waterp", "CHG", "CHGp", "Full", "Fullp"))

plot_richness(ps.Aim2, x = "Timepoint", measures=c("Shannon"), color = "WaterCHG") + 
  geom_boxplot(size = 1)  + 
  geom_point( size = 2.5)+
  facet_grid(rows = vars(Pig), col = vars(PmaxxGroup), scales = "free_x") + 
  scale_color_manual(values =c("#009EA3", "#F8CE62","#FD7A5E"))+
  theme_bw()

plot_richness(ps.Aim2, x = "Timepoint", measures=c("Shannon"), color = "WaterCHG", shape = "Pig") + 
  geom_boxplot(size = 1)  + 
  geom_point( size = 2.5)+
  facet_grid(col = vars(PmaxxGroup), scales = "free_x") + 
  scale_color_manual(values =c("#009EA3", "#F8CE62","#FD7A5E"))+
  theme_bw()

tab <-microbiome::alpha(ps.Aim2, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                   "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                   "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
#write.csv(tab, "./Aim2alphaTable.csv")


# Figure 2: relative abundnce tables 
OtherAim2rel <- as.data.frame(psmelt(ps.Aim2_REl))
OtherAim2 <- OtherAim2rel
OtherAim2 <- subset.data.frame(OtherAim2, Abundance >= 0.0001)
OtherAim2 <- OtherAim2%>% mutate(Genus = ifelse(Abundance < 0.02 , "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherAim2 <- OtherAim2 %>% mutate(Phylum = ifelse(Abundance < 0.02 , "Other", Phylum)) 
OtherAim2[OtherAim2 == "[Ruminococcus]"] <- "Ruminococcus"
OtherAim2[OtherAim2 == "Propionibacterium"] <- "Cutibacterium"
OtherAim2[OtherAim2 == "Post Tx"] <- "0 hr"
OtherAim2$Timepoint <- factor(OtherAim2$Timepoint, levels = c("Baseline", "0 hr","6 hr","12 hr", "24 hr", "48 hr"))
OtherAim2$Base.WaterCHG <- factor(OtherAim2$Base.WaterCHG, levels = c("Baseline", "Water","CHG","Full"))
table(OtherAim2$Genus, OtherAim2$Phylum)

OtherAim2$Genus<- factor(OtherAim2$Genus, levels = c("Other",
                                                     "Actinomadura","Actinomyces", "Arthrobacter", "Bifidobacterium", "Corynebacterium","Cutibacterium","Geodermatophilus", "Kocuria","Microbacterium", 
                                                     "Micrococcus", "Modestobacter", "Pseudonocardia", "Rothia","Streptomyces", #Actinobacteria / actinomycetota (14)
                                                     
                                                     "Aerococcus","Anaerococcus","Anaerovibrio", "Bacillus", "Clostridium","Coprococcus","Dorea", "Enterococcus","Epulopiscium","Faecalibacterium","Finegoldia", "Geobacillus",
                                                     "Granulicatella","Guggenheimella", "Jeotgalicoccus","Lactobacillus", "Lactococcus","Leuconostoc","Macrococcus","Oscillospira","Phascolarctobacterium", "Ruminococcus","Rummeliibacillus",
                                                     "Saccharibacillus", "Sarcina","SMB53","Sporosarcina", "Staphylococcus", "Streptococcus","Thermicanus","Trichococcus", "Turicibacter","Veillonella",  # Firmicutes/ Bacillota (31)
                                                     
                                                     "Flavobacterium","Hymenobacter","Parabacteroides", "Pedobacter","Prevotella","Sphingobacterium","Spirosoma", # Bacteroidetes / Bacteroidota (7)
                                                     "Acholeplasma",  #Mycoplasmatota/ Tenericutes
                                                     
                                                     "Acinetobacter","Actinobacillus","Alishewanella", "Brevundimonas","Delftia","Devosia","Diaphorobacter", "Enhydrobacter","Enterobacter","Erwinia",
                                                     "Janthinobacterium","Lysobacter","Massilia", "Moraxella","Proteus","Pseudomonas", "Psychrobacter","Rhodobacter", "Rubellimicrobium","Stenotrophomonas", # Proteobacteria /pseudomonadota (20)
                                                     
                                                     "Akkermansia")) #Verrucomicrobiota


ggplot(OtherAim2, aes(x=Timepoint, y = Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_grid(row=vars(PMAxx), col = vars(Base.WaterCHG), scales = "free_x")+
  ylab("Relative Abundance")+
  theme_bw()+
  scale_fill_manual(values = c("#EEEEEE",
                               colorRampPalette(c("#FFF1D7","#FED486","#FDB632","#F99536","#F47439","#D55E28","#9A3709" ))(12),
                               colorRampPalette(c("#DDFAF9","#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(26),
                               colorRampPalette(c("#C3ECFD","#64B3D9","#337799","#023B58","#022B40" ,"#011B28"))(17)))+
  theme(axis.text.x = element_text(size = 7)) +
  theme_light()+
  ggtitle("Genus present > 0.2% of reads in a sample")



# Figure 3B: Dot plot for key taxa in the Water CHG -viable communities
gen <- tax_glom(ps.Aim2, "Genus")
relAll.TV <- microbiome::transform(gen, "compositional")
relAll.V <- subset_samples(relAll.TV, PMAxx == "PMAxx")
relAll.V@sam_data$WaterCHG.Timepoint <- paste(relAll.V@sam_data$WaterCHG, relAll.V@sam_data$Timepoint) 
relAll.V.Grouped <- merge_samples2(relAll.V, "WaterCHG.Timepoint",  fun_otu = mean)
Aim2.Grouped <- psmelt(relAll.V.Grouped)
Aim2.Grouped  <- as.data.frame(Aim2.Grouped)
table(Aim2.Grouped$Phylum)
Aim2.Grouped$WaterCHG <-factor(Aim2.Grouped$WaterCHG, levels = c("Water", "CHG", "Full"))
Aim2.Grouped$Phylym<-factor(Aim2.Grouped$Phylum, levels = c("Other", "Actinobacteria", "Firmicutes", "Bacteroidota", "Proteobacteria"))
Aim2.Grouped$Timepoint <- factor(Aim2.Grouped$Timepoint, levels = c("Baseline", "0 hr","6 hr", "12 hr","24 hr", "48 hr"))

Aim2.Grouped.AllTop <- subset.data.frame(Aim2.Grouped, Genus %in% c("Acinetobacter", "Aerococcus", "Bacillus", "Corynebacterium", 
                                                                    "Macrococcus", "Moraxella", "Lactobacillus","Proteus","Pseudomonas",
                                                                    "Rothia", "SMB53","Staphylococcus", "Streptococcus", "Streptomyces", "Turicibacter"))

ggplot(Aim2.Grouped.AllTop, aes(x= Timepoint, y =Genus, Phylum, fill = Phylum, color = Phylum))+
  geom_point(aes(size = Abundance, fill = Phylum, color = Phylum)) +
  facet_grid(cols = vars(WaterCHG), scales = "free")+
  ylab("Mean relative Abundance")+
  scale_size(range=c(-1,9), breaks=c(0,0.01,0.05, .1,.2,.3, 0.4,0.5,0.6))+
  scale_x_discrete(labels = c("Baseline", "0 hr","6 hr", "12 hr","24 hr", "48 hr"))+
  scale_y_discrete(limits = rev(c("Corynebacterium", "Rothia" , "Streptomyces",
                                  "Aerococcus", "Bacillus", "Lactobacillus", "Macrococcus", "SMB53", "Staphylococcus", "Streptococcus", "Turicibacter",
                                  "Acinetobacter", "Moraxella", "Proteus", "Pseudomonas")))+
  theme_light()+
  scale_fill_manual(values = c("#FDB632", "#027979", "#023B58"))+ 
  scale_color_manual(values =c("#FDB632", "#027979", "#023B58"))+ #"#F47439"
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, vjust = 1), axis.text.y = element_text(size = 8)) +
  ggtitle("Mean Realtive abundance of genera")



# Supplemental table 4:beta diversity all (total and viable samples)
ps.Aim2.GenusGlom <- tax_glom(ps.Aim2, "Genus")
min(sample_sums(ps.Aim2.GenusGlom))# minimum sample read is 0
median(sample_sums(ps.Aim2.GenusGlom)) # 14809
TabS<- table(ps.Aim2.GenusGlom@sam_data$PMAxxSample, sample_sums(ps.Aim2.GenusGlom))

rarecurve(t(otu_table(ps.Aim2.GenusGlom)), step=100, ylim =c(0,60), xlim=c(0,3000))  ## considering using a read cut off of 5000 for beta diversity metrics 
ps.rareAllPigs = rarefy_even_depth(ps.Aim2.GenusGlom, rngseed=1, sample.size=1500, replace=F) # THIS IS THE ONE TO GO WITH for P1 

sampledf <- data.frame(sample_data(ps.rareAllPigs) )
Aim2_bray <- phyloseq::distance(ps.rareAllPigs, method = "bray")
adonis2(Aim2_bray ~ Pig, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 


### Supplemental Figure 3A-D and Supplemental table 5 - Pmaxx assessment 
GP = ps.rareAllPigs
GP.ord <- ordinate(GP, "NMDS", "bray")

Base.rare <-subset_samples(ps.rareAllPigs, Timepoint == "Baseline")
GP = Base.rare
GP.ord <- ordinate(GP, "NMDS", "bray")
plot_ordination(GP, GP.ord, type="samples", color="PMAxx", shape = "Pig") + 
  geom_point(size=6, alpha = 0.75) + ggtitle("Bray Curtis for all pigs Baseline Post treatemtnt timepoints ")+
  scale_color_manual(values = c("#023B58","#E3AC04"))+
  theme_bw()
sampledf <- data.frame(sample_data(Base.rare) )
Aim2_bray <- phyloseq::distance(Base.rare, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

PostTX.rare<-subset_samples(ps.rareAllPigs, Timepoint != "Baseline")
Water.048 <- subset_samples(PostTX.rare, WaterCHG == "Water")
GP = Water.048
GP.ord <- ordinate(GP, "NMDS", "bray")
plot_ordination(GP, GP.ord, type="samples", color="PMAxx", shape = "Pig") + 
  geom_point(size=6, alpha = 0.75) + ggtitle("Bray Curtis for all pigs Water Post treatemtnt timepoints ")+
  scale_color_manual(values = c("#023B58","#E3AC04"))+
  theme_bw()
sampledf <- data.frame(sample_data(Water.048) )
Aim2_bray <- phyloseq::distance(Water.048, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

CHG.048 <- subset_samples(PostTX.rare, WaterCHG == "CHG")
GP = CHG.048
GP.ord <- ordinate(GP, "NMDS", "bray")
plot_ordination(GP, GP.ord, type="samples", color="PMAxx", shape = "Pig") + 
  geom_point(size=6, alpha = 0.75) + ggtitle("Bray Curtis for all pigs CHG Post treatemtnt timepoints ")+
  scale_color_manual(values = c("#023B58","#E3AC04"))+
  theme_bw()
sampledf <- data.frame(sample_data(CHG.048) )
Aim2_bray <- phyloseq::distance(CHG.048, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

Full.048 <- subset_samples(PostTX.rare, WaterCHG == "Full")
GP = Full.048
GP.ord <- ordinate(GP, "NMDS", "bray")
plot_ordination(GP, GP.ord, type="samples", color="PMAxx", shape = "Pig") + 
  geom_point(size=6, alpha = 0.75) + ggtitle("Bray Curtis for all pigs Full Post treatemtnt timepoints ")+
  scale_color_manual(values = c("#023B58","#E3AC04"))+
  theme_bw()
sampledf <- data.frame(sample_data(Full.048) )
Aim2_bray <- phyloseq::distance(Full.048, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 


Water.all <- subset_samples(ps.rareAllPigs, WaterCHG == "Water")
PA <-Water.all
Base<-subset_samples(PA, Timepoint == "Baseline")
sampledf <- data.frame(sample_data(Base) )
Aim2_bray <- phyloseq::distance(Base, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.0<-subset_samples(PA, Timepoint == "0 hr")
sampledf <- data.frame(sample_data(HR.0) )
Aim2_bray <- phyloseq::distance(HR.0, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.6<-subset_samples(PA, Timepoint == "6 hr")
sampledf <- data.frame(sample_data(HR.6) )
Aim2_bray <- phyloseq::distance(HR.6, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.12<-subset_samples(PA, Timepoint == "12 hr")
sampledf <- data.frame(sample_data(HR.12) )
Aim2_bray <- phyloseq::distance(HR.12, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.24<-subset_samples(PA, Timepoint == "24 hr")
sampledf <- data.frame(sample_data(HR.24) )
Aim2_bray <- phyloseq::distance(HR.24, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.48<-subset_samples(PA, Timepoint == "48 hr")
sampledf <- data.frame(sample_data(HR.48) )
Aim2_bray <- phyloseq::distance(HR.48, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 


CHG.all <- subset_samples(ps.rareAllPigs, WaterCHG == "CHG")
PA <-CHG.all
Base<-subset_samples(PA, Timepoint == "Baseline")
sampledf <- data.frame(sample_data(Base) )
Aim2_bray <- phyloseq::distance(Base, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.0<-subset_samples(PA, Timepoint == "0 hr")
sampledf <- data.frame(sample_data(HR.0) )
Aim2_bray <- phyloseq::distance(HR.0, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.6<-subset_samples(PA, Timepoint == "6 hr")
sampledf <- data.frame(sample_data(HR.6) )
Aim2_bray <- phyloseq::distance(HR.6, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.12<-subset_samples(PA, Timepoint == "12 hr")
sampledf <- data.frame(sample_data(HR.12) )
Aim2_bray <- phyloseq::distance(HR.12, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.24<-subset_samples(PA, Timepoint == "24 hr")
sampledf <- data.frame(sample_data(HR.24) )
Aim2_bray <- phyloseq::distance(HR.24, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.48<-subset_samples(PA, Timepoint == "48 hr")
sampledf <- data.frame(sample_data(HR.48) )
Aim2_bray <- phyloseq::distance(HR.48, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 


Full.all <- subset_samples(ps.rareAllPigs, WaterCHG == "Full")
PA <-Full.all
Base<-subset_samples(PA, Timepoint == "Baseline")
sampledf <- data.frame(sample_data(Base) )
Aim2_bray <- phyloseq::distance(Base, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.0<-subset_samples(PA, Timepoint == "0 hr")
sampledf <- data.frame(sample_data(HR.0) )
Aim2_bray <- phyloseq::distance(HR.0, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.6<-subset_samples(PA, Timepoint == "6 hr")
sampledf <- data.frame(sample_data(HR.6) )
Aim2_bray <- phyloseq::distance(HR.6, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.12<-subset_samples(PA, Timepoint == "12 hr")
sampledf <- data.frame(sample_data(HR.12) )
Aim2_bray <- phyloseq::distance(HR.12, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.24<-subset_samples(PA, Timepoint == "24 hr")
sampledf <- data.frame(sample_data(HR.24) )
Aim2_bray <- phyloseq::distance(HR.24, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

HR.48<-subset_samples(PA, Timepoint == "48 hr")
sampledf <- data.frame(sample_data(HR.48) )
Aim2_bray <- phyloseq::distance(HR.48, method = "bray")
adonis2(Aim2_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+PMAxx, by= "margin", data = sampledf, permutations = 9999) # 


## Figure 3A - Viable microbiome only - ALL PIGS (with baseline indicated separately)
ps.Aim2.GenusGlom@sam_data$Timepoint <- factor(ps.Aim2.GenusGlom@sam_data$Timepoint , levels = c("Baseline","0 hr","6 hr","12 hr", "24 hr", "48 hr" ))
ps.Aim2.GenusGlom@sam_data$WaterCHG <- factor(ps.Aim2.GenusGlom@sam_data$WaterCHG , levels = c("Water", "CHG", "Full"))
ps.Aim2.GenusGlom@sam_data$Base.WaterCHG <- factor(ps.Aim2.GenusGlom@sam_data$Base.WaterCHG , levels = c("Baseline", "Water", "CHG", "Full"))

ps.Viable <- subset_samples(ps.Aim2.GenusGlom, PMAxx == "PMAxx")
table(sample_sums(ps.Viable))
rarecurve(t(otu_table(ps.Viable)), step=100, ylim =c(0,60), xlim=c(0,3000))  ## considering using a read cut off of 5000 for beta diversity metrics 
rare.Viable.AllP = rarefy_even_depth(ps.Viable , rngseed=1, sample.size=1500, replace=F) # THIS IS THE ONE TO GO WITH for P1 

GP =rare.Viable.AllP 
GP.ord <- ordinate(GP, "NMDS", "bray")
plot_ordination(GP, GP.ord, type="samples", color="Base.WaterCHG", shape = "Pig") + 
  geom_point(size=6, alpha = 0.75) + ggtitle("Bray Curtis for All pigs- V and T")+
  scale_color_manual(values = c("#CCCCCC", "#027979","#FDB632","#F47439"))+
  theme_light()


## Supplemental Figure 6A and supplemental table 7: distance from centroids for the Water- CHG compairison 
GP =rare.Viable.AllP 
BWCF.meta =as.data.frame(sample_data(rare.Viable.AllP))
BWCF <- factor( rare.Viable.AllP@sam_data$Base.WaterCHG)
V.dist <-distance(rare.Viable.AllP, method="bray")
dist.cent <- betadisper(V.dist, BWCF, type = "centroid", bias.adjust = F, sqrt.dist = F)
dist.cent
plot(dist.cent , main="PCoA")
boxplot(dist.cent, main="Distance to centroids")
TukeyHSD(dist.cent, ordered =T,
         conf.level = 0.95)

BWCF.dist<-as.data.frame(dist.cent$distances)
BWCF.dist2<-merge(BWCF.dist,BWCF.meta, by = 'row.names' )
BWCF.dist2$Base.WaterCHG<- factor(BWCF.dist2$Base.WaterCHG, levels =c("Baseline", "Water", "CHG", "Full" ))
BWCF.sum<- BWCF.dist2 %>% 
  group_by(Base.WaterCHG) %>% 
  summarise( mean = mean(`dist.cent$distances`),std = sd(`dist.cent$distances`))


ggplot(BWCF.dist2, aes(x= Base.WaterCHG, y= as.numeric(`dist.cent$distances`), color = Base.WaterCHG, fill = Base.WaterCHG ))+
  geom_boxplot(alpha = 0.7, linewidth = 2)+
  scale_y_continuous(limits = c(0,0.8))+
  scale_color_manual(values = c("#A7A9AC", "#012828","#B67702","#9A3709"))+
  scale_fill_manual(values = c("#EEEEEE", "#027979","#FDB632","#F47439"))+
  theme_light()


#Supplemental Table 6:  Viable Microbiome - compairisons of treatment groups at each timepoint 
PA <-rare.Viable.AllP 
sampledf <- data.frame(sample_data(PA) ) # all timepoints
Aim2_bray <- phyloseq::distance(PA, method = "bray")
adonis2(Aim2_bray ~ Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) 

Base<-subset_samples(PA, Timepoint == "Baseline")
sampledf <- data.frame(sample_data(Base) )
Aim2_bray <- phyloseq::distance(Base, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.0.48<-subset_samples(PA, Timepoint != "Baseline")
sampledf <- data.frame(sample_data(HR.0.48) )
Aim2_bray <- phyloseq::distance(HR.0.48, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.0<-subset_samples(PA, Timepoint == "0 hr")
sampledf <- data.frame(sample_data(HR.0) )
Aim2_bray <- phyloseq::distance(HR.0, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.6<-subset_samples(PA, Timepoint == "6 hr")
sampledf <- data.frame(sample_data(HR.6) )
Aim2_bray <- phyloseq::distance(HR.6, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.12<-subset_samples(PA, Timepoint == "12 hr")
sampledf <- data.frame(sample_data(HR.12) )
Aim2_bray <- phyloseq::distance(HR.12, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.24<-subset_samples(PA, Timepoint == "24 hr")
sampledf <- data.frame(sample_data(HR.24) )
Aim2_bray <- phyloseq::distance(HR.24, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.48<-subset_samples(PA, Timepoint == "48 hr")
sampledf <- data.frame(sample_data(HR.48) )
Aim2_bray <- phyloseq::distance(HR.48, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 


# Total Microbiome (not shown in paper)- compairisons of treatment groups at each timepoint 
ps.Total <- subset_samples(ps.Aim2.GenusGlom, PMAxx != "PMAxx")
table(sample_sums(ps.Total))
rarecurve(t(otu_table(ps.Total)), step=100, ylim =c(0,60), xlim=c(0,3000))  ## considering using a read cut off of 5000 for beta diversity metrics 
rare.Total.AllP = rarefy_even_depth(ps.Total , rngseed=1, sample.size=2000, replace=F) # THIS IS THE ONE TO GO WITH for P1 

PA <-rare.Total.AllP 
sampledf <- data.frame(sample_data(PA) ) # all timepoints
Aim2_bray <- phyloseq::distance(PA, method = "bray")
adonis2(Aim2_bray ~ Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+Base.WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

Base<-subset_samples(PA, Timepoint == "Baseline")
sampledf <- data.frame(sample_data(Base) )
Aim2_bray <- phyloseq::distance(Base, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.0<-subset_samples(PA, Timepoint == "0 hr")
sampledf <- data.frame(sample_data(HR.0) )
Aim2_bray <- phyloseq::distance(HR.0, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.6<-subset_samples(PA, Timepoint == "6 hr")
sampledf <- data.frame(sample_data(HR.6) )
Aim2_bray <- phyloseq::distance(HR.6, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.12<-subset_samples(PA, Timepoint == "12 hr")
sampledf <- data.frame(sample_data(HR.12) )
Aim2_bray <- phyloseq::distance(HR.12, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.24<-subset_samples(PA, Timepoint == "24 hr")
sampledf <- data.frame(sample_data(HR.24) )
Aim2_bray <- phyloseq::distance(HR.24, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 

HR.48<-subset_samples(PA, Timepoint == "48 hr")
sampledf <- data.frame(sample_data(HR.48) )
Aim2_bray <- phyloseq::distance(HR.48, method = "bray")
adonis2(Aim2_bray ~ WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(Aim2_bray ~ Pig+WaterCHG, by= "margin", data = sampledf, permutations = 9999) # 





# MAASLIN ASSESSMENTS
input_genus<- read.csv("./AllPigs_16s_Genus_Pivot.csv")
input_genus_data <- as.data.frame(input_genus )
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <-Aim2.meta
input_meta_file <-subset.data.frame(input_meta_file, Pig!= "NegEx")
input_meta_file$Time<-as.numeric(input_meta_file$Time)


# Supplemental fig 3 - PMAxx v. none For CHG and Full treatment groups only at all post intervention timepoints (0hr to 48hr)
V.T.0.48<-subset.data.frame(input_meta_file, Timepoint != "Baseline")
V.T.CHG.0.48 <- subset.data.frame(V.T.0.48, WaterCHG %in% c("CHG", "Full"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.T.CHG.0.48, 
  output = "./MAASLIN/PMAxxVnone.AllPigs.CHG.Full.0hrto48hr", 
  min_abundance = 0.01,
  fixed_effects = c("PMAxx"),
  reference= c("PMAxx", "Untreated"),
  random_effects = c("Pig"))



#Supplemental Figure 4 Effect of enviornment (water treated group)
Viable.Pigs<- subset.data.frame(input_meta_file, PMAxx == "PMAxx")
Viable.Water <- subset.data.frame(Viable.Pigs, WaterCHG == "Water")
Viable.Water.P1 <-subset.data.frame(Viable.Water, Pig == "P1")
Viable.Water.P2<-subset.data.frame(Viable.Water, Pig == "P2")
Viable.Water.P3<- subset.data.frame(Viable.Water, Pig == "P3")
Viable.Water.B.0<-subset.data.frame(Viable.Water, Timepoint %in% c("Baseline", "0 hr"))
Viable.Water.0.48<-subset.data.frame(Viable.Water, Timepoint !="Baseline")

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.Water.B.0, 
  output = "./MAASLIN/WaterBaselintTo0hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference= c("Timepoint", "Baseline"),
  random_effects = c("Pig"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.Water.0.48, 
  output = ".MAASLIN/Water.PostTreatment.OverTime_Cor_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("Time"),
  random_effects = c("Pig"))


###Supplemental Figure 5 - CHG groups (local CHG and full CHG prep) over time 
Viable.CHG <- subset.data.frame(Viable.Pigs, WaterCHG == "CHG")
Viable.CHG.B.0<-subset.data.frame(Viable.CHG, Timepoint %in% c("Baseline", "0 hr"))
Viable.CHG.0.48 <-subset.data.frame(Viable.CHG, Timepoint != "Baseline")

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.CHG.B.0, 
  output = "./MAASLIN/CHGBaselintTo0hr_allPigs", 
  min_abundance = 0.001,
  fixed_effects = c("Timepoint"),
  reference= c("Timepoint", "Baseline"),
  random_effects = c("Pig"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.CHG.0.48, 
  output = "./MAASLIN/CHG.PostTreatment.OverTime_Cor_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("Time"),
  random_effects = c("Pig"))


Viable.Full <- subset.data.frame(Viable.Pigs, WaterCHG == "Full")
Viable.Full.B.0<-subset.data.frame(Viable.Full, Timepoint %in% c("Baseline", "0 hr"))
Viable.Full.0.48 <-subset.data.frame(Viable.Full, Timepoint != "Baseline")

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.Full.B.0, 
  output = "./MAASLIN/FullBaselintTo0hr_allPigs", 
  min_abundance = 0.001,
  fixed_effects = c("Timepoint"),
  reference= c("Timepoint", "Baseline"),
  random_effects = c("Pig"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.Full.0.48, 
  output = "./MAASLIN/Full.PostTreatment.OverTime_Cor_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("Time"),
  random_effects = c("Pig"))


# Figure 3 and Supplemental Figure 6 - Water V, Chloroprep V. Fuill ALL post treatment timepoints combined 
Viable.0.48<- subset.data.frame(Viable.Pigs, Timepoint != "Baseline")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.0.48, 
  output = "./MAASLIN/Water.V.CHG.PostTreatmentAllTimepoints_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.0.48, 
  output = "./MAASLIN/Full.V.Water.V.CHG.PostTreatmentAllTimepoints_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Full"),
  random_effects = c("Pig"))

# Supplemental Table 8 continued. MAASLIN for CHG treatment groups v. eachother continued 
Viable.0<- subset.data.frame(Viable.Pigs, Timepoint == "0 hr")
Viable.0.P1 <- subset.data.frame(Viable.0, Pig == "P1")
Viable.0.P2 <- subset.data.frame(Viable.0, Pig == "P2")
Viable.0.P3 <- subset.data.frame(Viable.0, Pig == "P3")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.0, 
  output = "./MAASLIN/Water.V.CHG_0hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))
Viable.06<- subset.data.frame(Viable.Pigs, Timepoint == "6 hr")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.06, 
  output = "./MAASLIN/Water.V.CHG_06hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))
Viable.12<- subset.data.frame(Viable.Pigs, Timepoint == "12 hr")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.12, 
  output = "./MAASLIN/Water.V.CHG_12hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))
Viable.24<- subset.data.frame(Viable.Pigs, Timepoint == "24 hr")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.24, 
  output = "./MAASLIN/Water.V.CHG_24hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))
Viable.48<- subset.data.frame(Viable.Pigs, Timepoint == "48 hr")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Viable.48, 
  output = "./MAASLIN/Water.V.CHG_48hr_allPigs", 
  min_abundance = 0.01,
  fixed_effects = c("WaterCHG"),
  reference= c("WaterCHG", "Water"),
  random_effects = c("Pig"))

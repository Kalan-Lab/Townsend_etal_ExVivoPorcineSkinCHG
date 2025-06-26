library(broom)
library(ggrepel)
library(ggplot2)
library(tidyverse)
library(rstatix) 
library(xlsx)
library(RColorBrewer)
library(ANCOMBC) 
library(reshape2)
library(PairedData)

setwd(https://github.com/Kalan-Lab/Townsend_etal_ExVivoPorcineSkinCHG)

# Analysis of the Qubit DNA concentrations for the Pig CHG experiment (Not shown in paper)
QPC<- as.data.frame(./CompiledQubitDNA_PigCHGEx.v2.csv"))
QPC[QPC$DNAconc..ng.µl. == "Too low "] <- 0
QPC[QPC$DNAconc..ng.µl. == "too low "] <- 0
QPC$PmaxxGroup <- factor(QPC$PmaxxGroup, levels = c("Water", "Waterp", "CHG", "CHGp", "Full", "Fullp", "NegEx"))


ggplot(QPC, aes(x=PmaxxGroup, y=as.numeric(DNAconc..ng.µl.), color = WaterCHG))+
  geom_point(position = position_dodge2(width = 0.25))+
  scale_color_manual(values = c("#CA8402","#C1450B", "#CCCCCC","#011B28"))+
  theme_light() +
  ggtitle("DNA concentrations by Qubit") +
  ylab("DNA (ng/µl)")


QPC2 <- subset.data.frame(QPC, PmaxxGroup != "NegEx")
ANOVout<- aov(DNAconc..ng.µl. ~  PMAxx+WaterCHG +Time+Pig, data = QPC2) #two way anaova
summary(ANOVout)

PairedPC<- as.data.frame(read.csv("./PigCHGQubitPaired.csv")) # dataframe with samples in order 
PMAxx.df <- subset.data.frame(PairedPC,  PMAxx == "PMAxx")
PMAxx.vec <- PMAxx.df$DNAconc..ng.µl.
Untreated.df <- subset.data.frame(PairedPC,  PMAxx == "Untreated")
Untreated.vec <- Untreated.df$DNAconc..ng.µl.
wilcox_out <- wilcox.test(Untreated.vec, PMAxx.vec, paired = TRUE)
wilcox_out # p-value = 4.481e-12

# Analysis of Qubit DNA concentrations from the optomization experiments (Not Shown in Paper)
QPO<- as.data.frame(read.csv("./CompiledQubitDNA_PigOpt.csv"))
ANOVout<- aov(DNAconc..ng.µl. ~  PMAxx+LHKC +PMAxxConcentration+Pig, data = QPO) #two way anaova
summary(ANOVout)

KC <- subset.data.frame(QPO, LHKC != "Live" )
ANOVout<- aov(DNAconc..ng.µl. ~  PMAxx+LHKC +PMAxxConcentration+Pig, data = KC)
summary(ANOVout)

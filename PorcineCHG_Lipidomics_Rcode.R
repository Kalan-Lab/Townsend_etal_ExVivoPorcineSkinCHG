# Ex vivo porcine CHG project Lipidomics - R Code for analysis and figure generation 
library(phyloseq)
library(ggplot2)
library(ggforce)
library(microbiome)
library(tidyverse)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(DEqMS)
library(vegan)
library(mixOmics)
library(utils)
library(EnhancedVolcano)

setwd("/Users/liztown/Documents/KalanLab/Papers/PigSkinCHG/Rcode_FINAL")

#combining dataframes to make master metafile 
# lipid classes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3646453/
# ceramides https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2724059/

Lip.MF <- read.csv("./LipidClassMeta.csv")
Lip.MF<-as.data.frame(Lip.MF)
Lip.Class <- read.csv("./LipidClasses.csv")
Lip.Class<-as.data.frame(Lip.Class)
Lip.ID <-read.csv("./CompoundIDlist.csv")
Lip.ID <- as.data.frame(Lip.ID)

Lipid.Full <- merge(Lip.MF, Lip.Class, by.x = "LipidABV", by.y = "LipidClass.ABV", all = TRUE  )
Lipid.Full2 <- merge(Lip.ID, Lipid.Full, by.x = "CompoundName", by.y = "CompoundName", all = TRUE  )


#write.csv(Lipid.Full2, "./LipidClasses_Aim2_MetaFile.csv")

Lipid.FullMF <- read.csv("./LipidClasses_Aim2_MetaFile_Final.csv")
Lipid.FullMF<-unique(as.data.frame(Lipid.FullMF)) # full metafile with lipid classes 

# Figure 4 A relative abundance plots 
Pig2.Rel.lipids<- as.data.frame(read.csv("./Townsend_P2_NegPos_TissueMasNorm_RELATIVEabundance.csv"))
Pig3.Rel.lipids<- as.data.frame(read.csv("./Townsend_P3_NegPos_TissueMasNorm_RELATIVEabundance.csv"))

BothPigs.REL <- merge(Pig2.Rel.lipids, Pig3.Rel.lipids, by = "CompoundID", all = T )
row.names(BothPigs.REL) <- BothPigs.REL$CompoundID
BothPigs.REL<- BothPigs.REL[,-1]
BothPigs.REL.2<- setNames(melt(t(BothPigs.REL)), c('Sample','CompoundID', 'RelativeAbundance'))

PigMetadata<- read.csv("./FullPorcineSkinLipidomicManfest.csv")
PigMetadata2<-subset(PigMetadata, select = -c(PosPhaseSampleID,NegPhaseSampleID))
Lipid.FullMF <- read.csv("./LipidClasses_Aim2_MetaFile_Final.csv")
Lipid.FullMF2<-subset(Lipid.FullMF, select = -c(Formula, Mass, RT))
Lipid.FullMF2<-unique(as.data.frame(Lipid.FullMF2)) # full metafile with lipid classes 

BothPigs.REL.3 <- merge(BothPigs.REL.2,PigMetadata2, by.x="Sample", by.y="SampmleID", all = F )
BothPigs.REL.Full <- merge(BothPigs.REL.3,Lipid.FullMF2, by="CompoundID",all = F )
BothPigs.REL.Full<- subset.data.frame(BothPigs.REL.Full, RelativeAbundance > 0)
BothPigs.REL.Full$smp <- "All Samples"
BothPigs.REL.Full$Timepoint<-factor(BothPigs.REL.Full$Timepoint, levels = c("Baseline", "PostIntervention", "EndPoint"))
BothPigs.REL.Full$LIPIDClass_CommonName <- factor(BothPigs.REL.Full$LIPIDClass_CommonName, 
                                                  levels = c("Triglycerides", "Ceramides",
                                                             " Acyl Sterol Glycosides", "Fatty acids","Fatty amides","Fatty Esters",  # fatty Acyls
                                                             "Diglycerides","Glycosyldiradylglycerols",   # Glycerolipids
                                                             "Ether Lipids","Lysophospholipids","Phosphatidylcholines", "Phosphatidylethanolamines","Phosphatidylglycerols","Phosphatidylinositols","Phosphatidylserines", # Glycerophospholipids
                                                             "Hexosylceramides","Sphingoid bases" ,"Sphingomyelins",# sphingolipids
                                                             "Cholesteryl ester","Secosteroids","Sterols", "Steryl esters" # Sterol lipids 
                                                  ))

ggplot(BothPigs.REL.Full, aes(x = smp, y = RelativeAbundance, color = LIPIDClass_CommonName, fill = LIPIDClass_CommonName))+
  geom_bar(aes(y = RelativeAbundance), stat = "identity", position = "fill")+
  facet_grid(col = vars(Timepoint), scales = "free_x")+
  scale_fill_manual(values = c("#9A3709" ,"#64B3D9",
                               colorRampPalette(c("#FED486","#FDB632","#F99536","#F47439","#D55E28"))(6),
                               colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(7),
                               colorRampPalette(c("#337799","#023B58","#022B40" ,"#011B28"))(7)))+
  scale_color_manual(values = c("#9A3709" ,"#64B3D9",
                                colorRampPalette(c("#FED486","#FDB632","#F99536","#F47439","#D55E28"))(6),
                                colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(7),
                                colorRampPalette(c("#337799","#023B58","#022B40" ,"#011B28"))(7)))+
  theme_light()+
  ggtitle("Relative abundance of lipid species accross all samples by Timepoint")

#write.csv(BothPigs.REL.Full,"./BothPigs.RELATIVEABUNDANCE.Lipids2.csv")

low_abundance2<-as.data.frame(read.csv("./BothPigs.lowRELATIVEABUNDANCE.Lipids.csv"))
low_abundance <- subset.data.frame(low_abundance2, !LIPIDClass_CommonName %in% c("Ceramides", "Triglycerides"))
low_abundance$Timepoint <- factor(low_abundance$Timepoint, levels = c("Baseline", "PostIntervention", "EndPoint"))
low_abundance$LIPIDClass_CommonName <- factor(low_abundance$LIPIDClass_CommonName, 
                                              levels = c(" Acyl Sterol Glycosides", "Fatty acids","Fatty amides","Fatty Esters",  # fatty Acyls
                                                         "Diglycerides","Glycosyldiradylglycerols", "Triglycerides",  # Glycerolipids
                                                         "Ether Lipids","Lysophospholipids","Phosphatidylcholines", "Phosphatidylethanolamines","Phosphatidylglycerols","Phosphatidylinositols","Phosphatidylserines", # Glycerophospholipids
                                                         "Ceramides","Hexosylceramides","Sphingoid bases" ,"Sphingomyelins",# sphingolipids
                                                         "Cholesteryl ester","Secosteroids","Sterols", "Steryl esters" # Sterol lipids 
                                              ))
ggplot(low_abundance, aes(x = smp, y = RelativeAbundance, color = LIPIDClass_CommonName, fill = LIPIDClass_CommonName))+
  geom_bar(aes(y = RelativeAbundance), stat = "identity", position = "stack")+
  facet_grid(col = vars(Timepoint), scales = "free_x")+
  scale_y_continuous(limits = c(0, 0.1))+
  scale_fill_manual(values = c(colorRampPalette(c("#FED486","#FDB632","#F99536","#F47439","#D55E28"))(6),
                               colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(7),
                               colorRampPalette(c("#337799","#023B58","#022B40" ,"#011B28"))(7)))+
  scale_color_manual(values = c(colorRampPalette(c("#FED486","#FDB632","#F99536","#F47439","#D55E28"))(6),
                                colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(7),
                                colorRampPalette(c("#337799","#023B58","#022B40" ,"#011B28"))(7)))+
  theme_light()+
  ggtitle("Relative abundance of lipid species accross all samples")



#normalizing Lipidomics data in prep for remaining lipidomics figures and tables
# pig 2 data normalization
Pig2.lipids<- as.data.frame(read.csv("./Townsend_P2_NegPos_TissueMassNorm_input.csv"))
PigLipidMatrix <- subset(Pig2.lipids,select = -c(Mass,RT,CompoundName, Formula, Phase))
PigLipidMatrix[PigLipidMatrix== 0] <- 1   # pseudo count 
dat.lip <- as.matrix(PigLipidMatrix[,-1]) 
rownames(dat.lip) = PigLipidMatrix$CompoundID

countMat<-PigLipidMatrix
rownames(countMat) = countMat$CompoundID
countMat<- countMat[,-1]
countMat[countMat<= 1] <- 0
countMat[countMat > 0] <- 1
rs  <-  rowSums(countMat[ , ])  

dat.log = log2(dat.lip) # log transform data
dat.log = na.omit(dat.log) #remove rows with NAs
dat.log = equalMedianNormalization(dat.log) # median centring the data 
dat.log
boxplot(dat.log,las=2, cex.axis=.5)

# Pig 3 data normalization 
Pig3.lipids<- as.data.frame(read.csv("./Townsend_P3_NegPos_TissueMassNorm_input.csv"))
Pig3LipidMatrix <- subset(Pig3.lipids,select = -c(Mass,RT,CompoundName, Formula, Phase))
Pig3LipidMatrix[Pig3LipidMatrix== 0] <- 1   # pseudo count 
dat.lip3 <- as.matrix(Pig3LipidMatrix[,-1]) 
rownames(dat.lip3) = Pig3LipidMatrix$CompoundID

count3Mat<-Pig3LipidMatrix
rownames(count3Mat) = count3Mat$CompoundID
count3Mat<- count3Mat[,-1]
count3Mat[count3Mat<= 1] <- 0
count3Mat[count3Mat > 0] <- 1
rs3  <-  rowSums(count3Mat[ , ])  

dat.log3 = log2(dat.lip3) # log transform data
dat.log3 = na.omit(dat.log3) #remove rows with NAs
dat.log3 = equalMedianNormalization(dat.log3) # median centring the data 
dat.log3
boxplot(dat.log3,las=2, cex.axis=.5)



# Figure 4B and supplemental Figure 9
Pig2.log <- t(dat.log)
Pig3.log <- t(dat.log3)
both.pig.log <- merge(Pig2.log, Pig3.log, by = "row.names", all =T)
both.pig.log <-subset.data.frame(both.pig.log, Row.names != "P3_Water_EndPoint.2")
PigMetadata<- as.data.frame(read.csv("./FullPorcineSkinLipidomicManfest.csv"))
PigMetadata<-subset.data.frame(PigMetadata, SampmleID != "P3_Water_EndPoint.2")
PigMetadata$PigTimepoint <- paste(PigMetadata$Pig, PigMetadata$Timepoint)

result.pca.multi <- pca(both.pig.log, ncomp = 3)   # run the method
plotIndiv(result.pca.multi,
          group = PigMetadata$PigTimepoint,
          comp = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          pch=c(3,1,2,3,1,2),
          col =c("#B9B9B9", "#027978", "#FDB632",
                 "#B9B9B8", "#027977", "#FDB631"),
          # xlim = c(-63,82), ylim = c(-40,65),
          title = 'Both Pig lipids PCA plot') 
plotIndiv(result.pca.multi,
          group = PigMetadata$Pig,
          comp = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          pch=c(1,2),
          col =c("#FDB632","#023B58"),
          # xlim = c(-63,82), ylim = c(-40,65),
          title = 'Both Pig lipids PCA plot') 

plotIndiv(result.pca.multi,
          group = PigMetadata$Treatment.Base.Time,
          comp = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          pch=c(3,2,1,2,1,2,1),
          col =c("#B9B9B9", "#FDB632","#FDB631",
                 "#F47438","#F47439","#027978","#027979"),
          # xlim = c(-63,82), ylim = c(-40,65),
          title = 'Both Pig lipids PCA plot')

# Supplemental Table 10 - PERMANOVAS on PCA 
PCAcord<-result.pca.multi$variates[["X"]]
names <- rownames(PCAcord)
PCAcord1<- cbind(names,PCAcord)
colnames(PCAcord1)[colnames(PCAcord1) == "names"] = "SampmleID"
PCAcord2<-merge(PCAcord1, PigMetadata, by = "SampmleID")
PCAcord2$PC1 <- as.numeric(as.character(PCAcord2$PC1))
PCAcord2$PC2 <- as.numeric(as.character(PCAcord2$PC2))
sampledf <- data.frame(PigMetadata) 
PCAcord_adonis <- as.data.frame(PCAcord1[,-1])
PCAcord_adonis<- subset(PCAcord_adonis, select = -c(PC3))
PCAcord_adonis$PC1 <- as.numeric(as.character(PCAcord_adonis$PC1))
PCAcord_adonis$PC2 <- as.numeric(as.character(PCAcord_adonis$PC2))

adonis2(PCAcord_adonis ~ Pig, data = sampledf, method = "eu", by = "margin",permutations = 9999 )
adonis2(PCAcord_adonis ~ Treatment, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Base.Tx, data = sampledf, method = "eu", by = "margin",permutations = 9999 )
adonis2(PCAcord_adonis ~ Timepoint, data = sampledf, method = "eu",by = "margin", permutations = 9999 )

adonis2(PCAcord_adonis ~ Pig + Treatment, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Pig + Base.Tx, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Pig + Timepoint, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Timepoint + Pig, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Timepoint + Base.Tx, data = sampledf, method = "eu",by = "margin", permutations = 9999 )
adonis2(PCAcord_adonis ~ Timepoint + Treatment, data = sampledf, method = "eu", by = "margin",permutations = 9999 )


# Figure 4 C-E and Supplemental Table 11: differential expression analysis - difference between each timepoint 
Lipid.FullMF <- read.csv("./LipidClasses_Aim2_MetaFile_Final.csv")
Lipid.FullMF<-unique(as.data.frame(Lipid.FullMF))

Pig2.log <- dat.log
Pig3.log <- dat.log3
both.pig.log <- merge(Pig2.log, Pig3.log, by = "row.names", all =T)
row.names(both.pig.log) = both.pig.log$Row.names
both.pig.log <- both.pig.log[,-1]
both.pig.log <-subset(both.pig.log, select =-c(P3_Water_EndPoint.2))

countBothMat<-both.pig.log
SmallestValue<-countBothMat[which.min(countBothMat)] # - 19.777 (note this is log normalized data )
SmallestValue
countBothMat[is.na(countBothMat)] <- -20
countBothMat[countBothMat > -20] <- 1
countBothMat[countBothMat == -20] <- 0
rsBoth  <-  rowSums(countBothMat[ , ])  

colOrder <- colnames(both.pig.log)
colOrder 
cond = as.factor(c("Endpoint","Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt",
                   "Endpoint","Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt",
                   "Endpoint","Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt",
                   "Endpoint","Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt",
                   "Endpoint","Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt",
                   "Endpoint","Endpoint",
                   "Baseline","Baseline","Baseline",
                   "PostInt","PostInt","PostInt"
                   ))
design = model.matrix(~0+cond) # generating the design matrix. 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
x <- c("PostInt-Baseline","Endpoint-Baseline",
       "Endpoint-PostInt")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(both.pig.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

psm.count.table = data.frame(count = rsBoth , row.names =  rownames(countBothMat))
fit3$count = psm.count.table[rownames(fit3$coefficients), "count"]
fit4.b = spectraCounteBayes(fit3)
head(fit4.b)

PostInt.Baseline = outputResult(fit4.b,coef_col = 1) 
Endpoint.Baseline = outputResult(fit4.b,coef_col = 2) 
Endpoint.PostInt = outputResult(fit4.b,coef_col = 3) 

linker <- as.data.frame(Lipid.FullMF) 
PostInt.Baseline.res <- merge(PostInt.Baseline, linker, by.x="row.names", by.y = "CompoundID")
Endpoint.Baseline.res <- merge(Endpoint.Baseline, linker, by.x="row.names", by.y = "CompoundID")
Endpoint.PostInt.res <- merge(Endpoint.PostInt, linker, by.x="row.names", by.y = "CompoundID")

EnhancedVolcano(PostInt.Baseline.res , x = "logFC", y = "adj.P.Val", lab =PostInt.Baseline.res $CleanName, # BH adjusted pvalue 
                pCutoff = 1e-2, FCcutoff = 2 ,
               # ylim = c(0,5),
                xlab= "log2 Fold change - PostInt/ Baseline",
                pointSize = 4.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = ' 0hr/ Baseline')
EnhancedVolcano(Endpoint.Baseline.res, x = "logFC", y = "adj.P.Val", lab =Endpoint.Baseline.res$CleanName, # BH adjusted pvalue 
                pCutoff = 1e-2, FCcutoff = 2 ,
               # ylim = c(0,5),
                xlab= "log2 Fold change - 48hr/ Baseline ",
                pointSize = 4.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Water 48/ Baseline')

EnhancedVolcano(Endpoint.PostInt.res, x = "logFC", y = "adj.P.Val", lab =Endpoint.PostInt.res$CleanName, # BH adjusted pvalue 
                pCutoff = 1e-2, FCcutoff = 2 ,
                ylim = c(0,5),
                xlab= "log2 Fold change -  48 hr/ 0 hr",
                pointSize = 4.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = ' 48 hr/ 0 hr')

#write.csv(PostInt.Baseline.res, "./BaselineVPostInt0_DifferentialExpressionResults.csv")
#write.csv(Endpoint.Baseline.res ,"./BaselineVendpoint48_DifferentialExpressionResults.csv")
#write.csv(Endpoint.PostInt.res, "./PostInt0VEndpoint48_DifferentialExpressionResults.csv")


# Figure 4F- Key lipid changes on ex vivo skin over time
### (key lipids from the differential abundance analysis above [those with absolute(logFC >2) and adjusted p-value <0.01])
KeyLipids <- unique(c("PG 18:2 18:2","FA 20:3","FAHFA 18:1 20:3","FA 20:2","Cer NDS d19:0 30:1",
                      "Cer NDS d23:0 28:0","Cer ADS d51:0","FA 22:4","Cer EOS d26:1 23:0 18:2","Cer EOS d68:3",
                      "Cer EOS d18:1 30:0 18:2","Cer NS d21:1 30:0","Cer NDS d20:0 34:2","Cer NS d22:1 29:0","Cer NDS d19:0 32:1",
                      "Cer AP t52:1","Cer EODS d52:0","Cer EODS d18:0 15:1 20:0","PG 18:1 18:2","FA 19:0",
                      "FA 22:2","Cer NDS d16:0 16:0","PG 18:2 18:2","PS 18:0 20:3","Cer NS d21:1 30:0","PG 18:1 18:2","PC 38:5",
                      "PC 16:1 16:1","PE 14:1 16:1","PE 15:1 16:1","PC 18:1 18:2","PE 16:1 22:1","PI 16:1 16:1","SM d40:1",
                      "PC 32:2","PI 16:0 16:1","PC 40:2","PC 14:1 16:1","PC 15:0 16:1","PE 14:0 16:1","LPE 16:1","PC 16:0 16:1",
                      "TG 13:1 16:1 16:1","TG 45:3","TG 43:1","PE 35:1","PC 36:4","PE 38:4","PE 16:1 16:1","PE 17:0 16:1",
                      "FA 20:6","SM d42:2","TG 47:3","LPE 16:0","PC 36:4","SM d38:1","PE 15:0 16:1","PI 18:1 18:1","FA 20:5","PS 16:1 18:3",
                      "PE 16:1 18:1","PC 18:2 20:4","TG 46:5","PE 16:0 16:1","SM d42:1","LPC 16:1/0:0","TG 55:6","TG 13:1 16:0 16:1",
                      "CE 24:1","TG 16:0 16:0 18:1","PC 17:0 16:1","PC 16:1 20:4","TG 15:0 16:1 16:1","TG 48:5","PE 18:1 18:1",
                      "TG 44:4","PI 18:1 20:4","PC 40:2","PS 18:0 20:3","PE 35:1","PC 38:5","PE 14:1 16:1","PE 15:1 16:1",
                      "SM d40:1","PI 16:1 16:1","PC 32:2","PC 16:1 16:1","PC 15:0 16:1","PC 36:4","PE 38:4","PC 17:0 16:1",
                      "PC 14:1 16:1","LPE 16:1","SM d38:1","LPE 16:0","FA 20:6","PC 16:0 16:1","DG 50:9","PI 16:0 16:1",
                      "TG 55:6","PE 14:0 16:1","SM d42:1","TG 45:3","TG 13:1 16:1 16:1"))

Key.Lip.d<-as.data.frame(read.csv("./BothPigs.Key_byTime_RELATIVEABUNDANCE.Lipids.csv"))
row.names(Key.Lip.d) <- Key.Lip.d$CleanName
Key.Lip.d<- Key.Lip.d[,-1]
Key.Lip.d<- setNames(melt(t(Key.Lip.d)), c('Timepoint','CleanName', 'MeanRelativeAbundance'))
Lipid.FullMF3<-subset(Lipid.FullMF2, select= -c(CompoundID, CompoundName))
Key.Lip.d <- merge(Key.Lip.d,Lipid.FullMF3, by="CleanName",all = F )
Key.Lip.abd<- subset.data.frame(Key.Lip.d, CleanName %in%KeyLipids)
Key.Lip.abd$Timepoint <- factor(Key.Lip.abd$Timepoint, levels = c("Baseline", "PostIntervention", "EndPoint"))

Key.Lip.abd$LIPIDClass_CommonName <- factor(Key.Lip.abd$LIPIDClass_CommonName, 
                                            levels = c(" Acyl Sterol Glycosides", "Fatty acids","Fatty amides","Fatty Esters",  # fatty Acyls
                                                       "Diglycerides","Glycosyldiradylglycerols", "Triglycerides",  # Glycerolipids
                                                       "Ether Lipids","Lysophospholipids","Phosphatidylcholines", "Phosphatidylethanolamines","Phosphatidylglycerols","Phosphatidylinositols","Phosphatidylserines", # Glycerophospholipids
                                                       "Ceramides","Hexosylceramides","Sphingoid bases" ,"Sphingomyelins",# sphingolipids
                                                       "Cholesteryl ester","Secosteroids","Sterols", "Steryl esters" # Sterol lipids 
                                            ))

ggplot(Key.Lip.abd, aes(x= Timepoint, y =CleanName, fill = LIPIDClass_CommonName, color = LIPIDClass_CommonName))+
  geom_point(aes(size = 100*MeanRelativeAbundance,  fill = LIPIDClass_CommonName, color = LIPIDClass_CommonName)) +
  facet_grid(row = vars(LargerLipidGroup.ABV), scales = "free_y", space = "free_y")+
  ylab("Mean log10 lipid Abundance")+
  scale_size(range=c(0,8), breaks = c(0,0.01,0.05, 0.1,0.25,0.5,0.75,1))+
  scale_y_discrete( limits = rev)+
  theme_light()+
  theme(axis.text.y = element_text(size = 4)) +
  scale_fill_manual(values = c("#FDB632","#C07E02","#F47439","#B14715",
                               colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(6),
                               "#3883A8","#1D5E7E", "#023954"))+
  scale_color_manual(values = c("#FDB632","#C07E02","#F47439","#B14715",
                                colorRampPalette(c("#88D2D2","#33AAAA","#027979","#024E51","#012228" ))(6),
                                "#3883A8","#1D5E7E","#023954"))+
  ggtitle("Mean % relative Key Lipid abundance at each timepoint")


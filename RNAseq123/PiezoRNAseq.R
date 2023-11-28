## RNAseq Analysis on Piezo RNA Seq using edgeR
## following: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(readr)
library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse)


##LOADING IN THE DATA##

path = "/Users/mtowriss/Desktop/Pz_RNA_seq/RNA-SEQ-FINAL"
setwd(path)
counts <- read.delim("R539_count.txt", row.names = 1)
head(counts)

#matrix gene names and ensembl IDs (column bind) 
genenames = cbind(counts$gene_name,rownames(counts))
genenames = as.data.frame(genenames)
colnames(genenames) = c("genesymbol","EnsemblID")

d0 <- DGEList(counts [,-1])
d0$genes = genenames

##getting rid of low expressed genes (cutoff for the number of samples that we want above 1, getting rid of things that are less than one group)
cutoff <- 3
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

#filtering plot
library(RColorBrewer)
nsamples <- ncol(d0)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)

#plot:
pdf('FilteringCPM_plots_cutoff3.pdf')
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(d0, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")

#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(d, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()

#Derive experiment information from the sample names
snames <- colnames(counts[,-1])
snames
sample<-substr(snames, 2, nchar(snames) - 8)
injection<-substr(snames, 3, nchar(snames) - 5)
sort<-substr(snames, 6, nchar(snames) - 0)
sample
injection
sort
#group interaction, set levels
group <- interaction(injection,sort)
group = factor(group,levels = c("PBS.Pzneg","PBS.Pzpos","LPS.Pzneg", "LPS.Pzpos"))
group

#make data frame for metadata
metadata = cbind(snames,group,sample,injection,sort)
metadata = data.frame(metadata)


#####CREATING THE MDS PLOTS#####
##creating a pdf for the MDS plots
#plot:
pdf('MDSplots.pdf')
par(mfrow=c(2,2))
#Multidimensional scaling (MDS) plot * colour by group
plotMDS(d, col = as.numeric(group))
title("SortxInjection")
#Multidimensional scaling (MDS) plot * colour by sample
plotMDS(d, col = as.numeric(sample))
title("Sample")
#Multidimensional scaling (MDS) plot * colour by injection
injection=factor(injection)
plotMDS(d, col = as.numeric(injection))
title("Injection")
#Multidimensional scaling (MDS) plot * colour by sort
sort=factor(sort)
plotMDS(d, col = as.numeric(sort))
title("Sort")

dev.off()

###Plot library size###
#plot library sizes
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(d$samples$lib.size,names=colnames(d),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

###Normalize Log counts###
# Get log2 counts per million
logcounts <- cpm(d,log=TRUE)
# Check distributions of samples using boxplots
pdf('NonNormalizedLogCPM1.pdf',w=30,h=10)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

####Design model matrix####
mm <- model.matrix(~0 + group)
##Run Calculation Normalization
DGE=calcNormFactors(d,method =c("TMM")) 
pdf('Voom.pdf',w=6,h=4)
v=voom(DGE,mm,plot=T)
dev.off()

#Duplicate correction because two samples from one animal condition
sample=factor(sample)
corfit <- duplicateCorrelation(v, mm, block=sample)
fit <- lmFit(v, mm, block=sample, correlation = corfit$consensus)


#####Average CPM for each condition#####
###get log2CPM counts from voom and put in dataframe:
library(plotrix)
#average log2 CPM and sem
countdf <- as.data.frame(v$E)
countdf$GeneID <- rownames(v$E)
#add gene names
DF <- merge(countdf,genenames, by.x ="GeneID",by.y="EnsemblID")
#write as csv
write.csv(DF,file="log2CPMvalues.csv")
#summarize 
countdf2 <- DF %>% group_by(GeneID,genesymbol) %>% tidyr::gather(sample,log2CPM, 2:(ncol(DF)-1)) 
countdf2 <- as.data.frame(countdf2)
metadata$magic=paste(metadata$injection,metadata$sort)
countdf3 <-merge(countdf2,metadata,by.y="snames",by.x="sample",all=T)

#Write CSV for Gene Summary
GeneSummary <- countdf3 %>% group_by(GeneID,genesymbol, magic) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  tidyr::spread(magic ,meanlog2CPM)
write.csv(GeneSummary, file = "AverageLog2CPM.csv")


###make contrastmatrix
contr.matrix <- makeContrasts(
  LPSPzposvsLPZPzneg = groupLPS.Pzpos - groupLPS.Pzneg,
  PBSPzposvsPBSPzneg = groupPBS.Pzpos - groupPBS.Pzneg,
  LPSPzposvsPBSPzpos = groupLPS.Pzpos - groupPBS.Pzpos,
  LPSPznegvsPBSPzneg = groupLPS.Pzneg - groupPBS.Pzneg,
  LPSPzposvsPBSPzneg = groupLPS.Pzpos - groupPBS.Pzneg,
  levels = colnames(mm))
contr.matrix

tmp <- contrasts.fit(fit, contr.matrix)
tmp <- eBayes(tmp)

pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(tmp, main="Final model: Mean variance trend")
dev.off()

#####MAKE CONTRASTS + GET DEGS#####
library(calibrate)

####Make contrasts between PBS.PZneg and PBS.Pzpos (upregulated in PBS Pz-positive)
PBSPzposvsPBSPzneg <- makeContrasts(groupPBS.Pzpos - groupPBS.Pzneg, levels = colnames(coef(fit)))
PBSPzposvsPBSPzneg.contr <- contrasts.fit(fit, PBSPzposvsPBSPzneg)
PBSPzposvsPBSPzneg.contr.ebayes <- eBayes(PBSPzposvsPBSPzneg.contr)
  #toptable
  top.table.PBSPzposvsPBSPzneg <- topTable(PBSPzposvsPBSPzneg.contr.ebayes, sort.by = "P", n = Inf)
  head(top.table.PBSPzposvsPBSPzneg, 20)
  #how many DEGS?
  length(which(top.table.PBSPzposvsPBSPzneg$adj.P.Val < 0.05))
  #savetop table as a file
  top.table.PBSPzposvsPBSPzneg$Gene <- rownames(top.table.PBSPzposvsPBSPzneg)
  write.table(top.table.PBSPzposvsPBSPzneg, file = "PBSPz+vsPBSPz-.txt", row.names = F, sep = "\t", quote = F)
##Make Volcano Plots for the contrast
  #volcano plot
  pdf(file = "PBSPz-vsPz+_Volcano.pdf", wi = 9, he = 6, useDingbats=F)
  with(top.table.PBSPzposvsPBSPzneg, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main="PBS Pz+ vs PBS Pz-", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(top.table.PBSPzposvsPBSPzneg, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(top.table.PBSPzposvsPBSPzneg, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  with(subset(top.table.PBSPzposvsPBSPzneg, adj.P.Val<0.05 & abs(logFC)>5), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.6))
  dev.off()

####Make contrasts between LPS.PZneg and LPS.Pzpos (upregulated in LPS Pz-positive)
LPSPzposvsLPZPzneg <- makeContrasts(groupLPS.Pzpos - groupLPS.Pzneg, levels = colnames(coef(fit)))
LPSPzposvsLPZPzneg.contr <- contrasts.fit(fit, LPSPzposvsLPZPzneg)
LPSPzposvsLPZPzneg.contr.ebayes <- eBayes(LPSPzposvsLPZPzneg.contr)
  #toptable
  top.table.LPSPzposvsLPSPzneg <- topTable(LPSPzposvsLPZPzneg.contr.ebayes, sort.by = "P", n = Inf)
  head(top.table.LPSPzposvsLPSPzneg, 20)
  #how many DEGS?
  length(which(top.table.LPSPzposvsLPSPzneg$adj.P.Val < 0.05))
  #savetop table as a file
  top.table.LPSPzposvsLPSPzneg$Gene <- rownames(top.table.LPSPzposvsLPSPzneg)
  write.table(top.table.LPSPzposvsLPSPzneg, file = "LPSPz+vsLPSPz-.txt", row.names = F, sep = "\t", quote = F)
##Make Volcano Plots for the contrast
  #volcano plot
  pdf(file = "LPSPz-vsPz+_Volcano.pdf", wi = 9, he = 6, useDingbats=F)
  with(top.table.LPSPzposvsLPSPzneg, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main="LPS Pz+ vs LPS Pz-", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(top.table.LPSPzposvsLPSPzneg, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(top.table.LPSPzposvsLPSPzneg, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  with(subset(top.table.LPSPzposvsLPSPzneg, adj.P.Val<0.05 & abs(logFC)>5), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.6))
  dev.off()

####Make contrasts between LPSPzpos vs PBSPzpos (upregulated in LPS Pz-positive)
LPSPzposvsPBSPzpos <- makeContrasts(groupLPS.Pzpos - groupPBS.Pzpos, levels = colnames(coef(fit)))
LPSPzposvsPBSPzpos.contr <- contrasts.fit(fit, LPSPzposvsPBSPzpos)
LPSPzposvsPBSPzpos.contr.ebayes <- eBayes(LPSPzposvsPBSPzpos.contr)
  #toptable
  top.table.LPSPzposvsPBSPzpos <- topTable(LPSPzposvsPBSPzpos.contr.ebayes, sort.by = "P", n = Inf)
  head(top.table.LPSPzposvsPBSPzpos, 20)
  #how many DEGS?
  length(which(top.table.LPSPzposvsPBSPzpos$adj.P.Val < 0.05))
  #savetop table as a file
  top.table.LPSPzposvsPBSPzpos$Gene <- rownames(top.table.LPSPzposvsPBSPzpos)
  write.table(top.table.LPSPzposvsPBSPzpos, file = "LPSPz+vsPBSPz+.txt", row.names = F, sep = "\t", quote = F)
##Make Volcano Plots for the contrast
  #volcano plot
  pdf(file = "LPSPz+vsPBSPz+_Volcano.pdf", wi = 9, he = 6, useDingbats=F)
  with(top.table.LPSPzposvsPBSPzpos, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main="LPS Pz+ vs PBS Pz+", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(top.table.LPSPzposvsPBSPzpos, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(top.table.LPSPzposvsPBSPzpos, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  with(subset(top.table.LPSPzposvsPBSPzpos, adj.P.Val<0.05 & abs(logFC)>5), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.6))
  dev.off()

####Make contrasts between LPSPzneg vs PBSPzneg (upregulated in LPS Pz-negative)
LPSPznegvsPBSPzneg <- makeContrasts(groupLPS.Pzneg - groupPBS.Pzneg, levels = colnames(coef(fit)))
LPSPznegvsPBSPzneg.contr <- contrasts.fit(fit, LPSPznegvsPBSPzneg)
LPSPznegvsPBSPzneg.contr.ebayes <- eBayes(LPSPznegvsPBSPzneg.contr)
  #toptable
  top.table.LPSPznegvsPBSPzneg <- topTable(LPSPznegvsPBSPzneg.contr.ebayes, sort.by = "P", n = Inf)
  head(top.table.LPSPznegvsPBSPzneg, 20)
  #how many DEGS?
  length(which(top.table.LPSPznegvsPBSPzneg$adj.P.Val < 0.05))
  #savetop table as a file
  top.table.LPSPznegvsPBSPzneg$Gene <- rownames(top.table.LPSPznegvsPBSPzneg)
  write.table(top.table.LPSPznegvsPBSPzneg, file = "LPSPz-vsPBSPz-.txt", row.names = F, sep = "\t", quote = F)
##Make Volcano Plots for the contrast
  #volcano plot
  pdf(file = "LPSPz-vsPBSPz-_Volcano.pdf", wi = 9, he = 6, useDingbats=F)
  with(top.table.LPSPznegvsPBSPzneg, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main="LPS Pz- vs PBS Pz-", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(top.table.LPSPznegvsPBSPzneg, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(top.table.LPSPznegvsPBSPzneg, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  with(subset(top.table.LPSPznegvsPBSPzneg, adj.P.Val<0.05 & abs(logFC)>5), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.6))
  dev.off()

####Make contrasts between LPSPzpos vs PBSPzneg (upregulated in LPS Pz-positive)
LPSPzposvsPBSPzneg <- makeContrasts(groupLPS.Pzpos - groupPBS.Pzneg, levels = colnames(coef(fit)))
LPSPzposvsPBSPzneg.contr <- contrasts.fit(fit, LPSPzposvsPBSPzneg)
LPSPzposvsPBSPzneg.contr.ebayes <- eBayes(LPSPzposvsPBSPzneg.contr)
  #toptable
  top.table.LPSPzposvsPBSPzneg <- topTable(LPSPzposvsPBSPzneg.contr.ebayes, sort.by = "P", n = Inf)
  head(top.table.LPSPzposvsPBSPzneg, 20)
  #how many DEGS?
  length(which(top.table.LPSPzposvsPBSPzneg$adj.P.Val < 0.05))
  #savetop table as a file
  top.table.LPSPzposvsPBSPzneg$Gene <- rownames(top.table.LPSPzposvsPBSPzneg)
  write.table(top.table.LPSPzposvsPBSPzneg, file = "LPSPz+vsPBSPz-.txt", row.names = F, sep = "\t", quote = F)
##Make Volcano Plots for the contrast
  #volcano plot
  pdf(file = "LPSPz+vsPBSPz-_Volcano.pdf", wi = 9, he = 6, useDingbats=F)
  with(top.table.LPSPzposvsPBSPzneg, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main="LPS Pz+ vs PBS Pz-", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(top.table.LPSPzposvsPBSPzneg, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(top.table.LPSPzposvsPBSPzneg, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  with(subset(top.table.LPSPzposvsPBSPzneg, adj.P.Val<0.05 & abs(logFC)>5), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.6))
  dev.off()



###MakeGlimmaPlots

################################################################################################
#GLIMMA interactive plot building
################################################################################################
#BiocManager::install("Glimma", version = "3.8")
#http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
library(Glimma)
  
#writes out html file:
glMDSPlot(v, groups=group)

dt <- decideTests(tmp)
summary(dt)
write.csv(summary(dt), file = "DEGcounts.csv")

for (COEF in 1:5) {
glMDPlot(tmp, counts=v$E,transform=FALSE,anno=tmp$genes,
         coef=COEF, status=dt, main=colnames(tmp)[COEF],
        groups=group, folder="glimma_results", launch=FALSE, html = paste("MD-Plot",
        colnames(contr.matrix)[COEF]))}

###Gene Set Enrichment

###Go Term Analysis - Kegg Pathways
##enricher website??? R version? 

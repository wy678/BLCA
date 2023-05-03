
rm(list=ls())
library(WGCNA)
library(data.table)
library(stringr)
library(gplots)
library(Biobase)
allowWGCNAThreads()

setwd("xx\\xxxx\\xxxx\\xxxxx\\xxxxxx")

load("xx/xxxx/xxxx/xxxx/risk_data.rda")
load("xx/xxxx/xxxx/xxxxx/VAEN_GDSC.A.pred_TCGA.rda")

rownames(risk_data)<-substr(rownames(risk_data),1,15)
exp<-GDSC[which(GDSC$Sample %in% rownames(risk_data)),]
risk_data<-risk_data[exp$Sample,]
rownames(exp)<-exp$Sample
exp<-exp[,-c(1,2)]


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(exp, powerVector = powers, verbose = 5)

pdf("1.软阈值筛选.pdf",height = 6,width = 10)

par(mfrow = c(1,2));
cex1 = 0.9;


plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()





ifelse(is.na(sft$powerEstimate) == T ,softPower <- 6 ,softPower <- sft$powerEstimate)
adjacency = adjacency(exp, power = softPower)

TOM = TOMsimilarity(adjacency)

dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

pdf("2.检验选定的β值.pdf",height = 6,width = 10)

ADJ1_cor <- abs(WGCNA::cor( exp,use = "p" ))^softPower

k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
 

par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()




net = blockwiseModules(
  exp,
  power = softPower,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  #saveTOMs = TRUE,
  #saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)


pdf("3.基因聚类模块.pdf",height = 6,width = 10)

mergedColors = labels2colors(net$colors)
table(mergedColors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

pdf("0.样本聚类.pdf",height = 6,width = 10)

datExpr_tree<-hclust(dist(exp), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
dev.off()



design<-model.matrix(~0+risk_data[,4])
colnames(design)<-c("High","Low")

moduleColors <- labels2colors(net$colors)

MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0);
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(exp))


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf("4.模块和性状的关系.pdf")
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),####
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(exp, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exp)))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")





Cluster<- as.data.frame(design[,2])
names(Cluster)<-"Subtype"
geneTraitSignificance = as.data.frame(cor(exp, Cluster, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exp)))
names(geneTraitSignificance) = paste("GS.", names(Cluster), sep="")
names(GSPvalue) = paste("p.GS.", names(Cluster), sep="")

pdf("5.blue_Cluster.pdf",height = 6,width = 6)


module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "NES for Cluster score",
                   main = paste("Module membership vs. Risk\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()



module_gene_MMGS=cbind(MM=abs(geneModuleMembership[which(moduleColors=="blue"), column]),
                       GS=abs(geneTraitSignificance[which(moduleColors=="blue"),1]))
rownames(module_gene_MMGS)=rownames(geneModuleMembership)[which(moduleColors=="blue")]
hub_gene1=rownames(module_gene_MMGS[which(module_gene_MMGS[,1]>0.5&module_gene_MMGS[,2]>0.4),])

save(module_gene_MMGS,file="blue_MMGS.rda")
write.table(hub_gene1,"hub模块(MM0.5,GS0.4).txt",quote = F,col.names = F,row.names = F)

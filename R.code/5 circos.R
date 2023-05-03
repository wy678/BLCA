
rm(list=ls())

library(HGNChelper)

setwd("x\\xxxx\\xxxx\\xxxx\\xxxxxx")

load("xx/xxxx/xxxx/xxxx/xxxxx/BLCA_all_data_by_genes.rda")
load("xx/xxxx/xxxx/xxxx/risk_data.rda")

data<-BLCA_CN_genes[,9:416]
rownames(risk_data)<-substr(gsub("-","\\.",rownames(risk_data)),1,16)
colnames(data)<-substr(colnames(data),1,16)
a<-risk_data[colnames(data,1,16),]
a<-na.omit(a)
data<-data[,rownames(a)]
data<-rbind(data,a$score)

dat=matrix(NA,nrow=nrow(data)-1,ncol=2)
for(i in 1:(nrow(data)-1))
{
  score=as.numeric(data[23345,])
  gene=as.numeric(data[i,])
  dat[i,1]=rownames(data)[i]
  dat[i,2]=cor(score,gene,method = "spearman")
  print(i)
  
}
head(dat)
dat[,2]=as.numeric(dat[,2])
colnames(dat)=c("Gene","Cor")

BLCA_CN<-merge(BLCA_CN_genes[,1:8],dat,by.x="GeneSymbol",by.y="Gene")
save(BLCA_CN,file = "BLCA_CN.rda")

############################################
rm(list=ls())

library(circlize)

setwd("xx\\xxx\\xxxx\\xxxx\\xxxxx")
load("BLCA_CN.rda")
load("xx/xxxx/xxxx/xxxx/TCGA_BLCA_edgeR_signif.rda")

nrDEG_edgeR_signif<-nrDEG_edgeR_signif[order(nrDEG_edgeR_signif$PValue),]
BLCA_CN<-BLCA_CN[which((BLCA_CN$GeneSymbol) %in% rownames(nrDEG_edgeR_signif)[1:500]),]

BLCA_CN$Cor<-as.numeric(BLCA_CN$Cor)
BLCA_CN$Chr<-paste0("chr",BLCA_CN$Chr)
BLCA_CN$Chr[which(BLCA_CN$Chr == "chr23")]<-"chrX"

a<-nrDEG_edgeR_signif[BLCA_CN$GeneSymbol,]
low<-rownames(a[which(a$logFC<0),])
high<-rownames(a[which(a$logFC>0),])

BLCA_CN$Cor1<-rep(NA,nrow(BLCA_CN))
BLCA_CN$Cor2<-rep(NA,nrow(BLCA_CN))
BLCA_CN[which(BLCA_CN$GeneSymbol %in% high),10]<-BLCA_CN[which(BLCA_CN$GeneSymbol %in% high),9]
BLCA_CN[which(BLCA_CN$GeneSymbol %in% low),11]<-BLCA_CN[which(BLCA_CN$GeneSymbol %in% low),9]

b<-BLCA_CN[,c(2,3,4,9,10,11)]


library(RColorBrewer)
mycol<-colorRampPalette(brewer.pal(12,'Set3'))(23)

col_fun <- colorRamp2(
  c(-0.3, 0, 0.3), 
  c("#00AFFF", "white", "#FC0E00")
)

c<-data.frame(matrix(NA,nrow = 23,ncol = 2))
c[,1]<-paste0("chr",c(1:22,"X"))
c[,2]<-mycol
linecol<-c[match(b$Chr,c[,1]),2]

pdf("Rcircos.pdf")

circos.initializeWithIdeogram(plotType = NULL)
circos.track(
  ylim = c(0, 1), 
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1)
    circos.text(
      mean(xlim), mean(ylim), chr, cex = 0.75,
      col = "white", facing = "inside",
      niceFacing = TRUE
    )
  }, 
  track.height = 0.15, bg.border = NA,bg.col =mycol
)
circos.genomicTrack(
  b, numeric.column = c(5,6),track.height = 0.15,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, numeric.column = 1,cex=0.5,col = "red")
    circos.genomicPoints(region, value, numeric.column = 2,cex=0.5,col = "blue")
  }
)
circos.genomicTrack(
  b, 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region, value, ytop.column = 1,
      ybottom = 0,
      col = ifelse(value[[1]] > 0, "#ef8a62", "#67a9cf"),
      ...
    )
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey50")
  }
)
circos.genomicHeatmap(
  b[,1:4], col = col_fun, side = "inside", connection_height = 0.05,heatmap_height = 0.1,
  border = "white",line_col = linecol
)
circos.clear()
dev.off()

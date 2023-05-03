
############################
rm(list=ls())

library(GSVA)
library(Biobase)
library(clusterProfiler)
library(tidyverse)

setwd("xx\\xxxx\\xxxx\\xxxx")
load("xx/xxxx/xxxx/TCGA_BLCA.rda")

#############################################
distCor <- function(x) as.dist(1-cor(x))
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  if (scale=="row") z <- t(scale(t(x))) else z<-x
  if (scale=="col") z <- scale(x)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  hcl_row <- hclust(distCor(t(z)), method=method)
  hcl_col <- hclust(distCor(z), method=method)
  return(list(data=z, Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
}

###############################################

hall<-read.gmt("G:\\work\\pain\\diff\\GSEA\\h.all.v7.5.1.symbols.gmt")
hall<-separate(data = hall, col = term, into = c("Hallmark", "term"), sep = "HALLMARK_")
genelist<-split(hall$gene,hall$term)

exp<-exprs(TCGA_BLCA)
GSVA<-gsva(exp,genelist,method="ssgsea",abs.ranking=F)
TCGA_GSVA<-zClust(GSVA)$data

save(TCGA_GSVA,file = "TCGA_GSVA.rda")

#########################################
rm(list=ls())

library(GEOquery)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(Hmisc)
setwd("xx\\xxxx\\xxxx\\xxxx")

load("xx/xxxx/xxxx/xxxx/TCGA_GSVA.rda")
load("xx/xxxx/xxxx/TCGA_BLCA.rda")
load("xx/xxxx/xxxx/xxxx/risk_data.rda")

risk_data$group<-capitalize(risk_data$group)
pdata<-pData(TCGA_BLCA)
pdata<-pdata[,c("OS.time","OS","age_at_initial_pathologic_diagnosis","gender","ajcc_pathologic_tumor_stage")]
pdata<-pdata[rownames(risk_data),]
pdata$Subtype<-risk_data$group
colnames(pdata)<-c("times","status","Age","Gender","Stage","Subtype")
pdata[,5]<-as.character(pdata[,5])
pdata<-pdata[-which(is.na(pdata$Stage) == T),]
pdata$Age<-as.numeric(pdata$Age)
pdata<-na.omit(pdata)
pdata[which(pdata[,3] > median(pdata[,3])),3]<-"Old"
pdata[which(pdata[,3] != "Old"),3]<-"Young"
pdata<-pdata[,-c(1,2)]

Cluster11<-rownames(pdata)[which(pdata$Subtype == "High")]
Cluster22<-rownames(pdata)[which(pdata$Subtype == "Low")]

NES_Cluster1<-TCGA_GSVA[,Cluster11]
NES_Cluster2<-TCGA_GSVA[,Cluster22]


Sigmark<-function(plot1){
  plot1[which(plot1[,1]>0.05|plot1[,1]==1),2]="ns"
  plot1[which(0.01<plot1[,1]&plot1[,1]<0.05),2]="*"
  plot1[which(0.001<plot1[,1]&plot1[,1]<0.01),2]="**"
  plot1[which(plot1[,1]<0.001),2]="***"
  return(plot1[,2])
}

random_wilcox<-data.frame()
for(i in 1:nrow(TCGA_GSVA)){
  random_wilcox[i,1]<-wilcox.test(NES_Cluster1[i,],NES_Cluster2[i,],
                                  exact=FALSE,correct=FALSE)$p.value
  random_wilcox[i,2]<-mean(NES_Cluster1[i,])
  random_wilcox[i,3]<-mean(NES_Cluster2[i,])
  if (random_wilcox[i,2]>random_wilcox[i,3]){
    random_wilcox[i,4]<-"High"
  }else if (random_wilcox[i,2]<=random_wilcox[i,3]){
    random_wilcox[i,4]<-"Low"
  }
}
random_wilcox[,5]<-Sigmark(random_wilcox)
row.names(random_wilcox)<-row.names(TCGA_GSVA)
colnames(random_wilcox)<-c("P-value","Cluster 1","Cluster 2","label","Signifi")

random_wilcox<-random_wilcox[order(random_wilcox$label),]

random_wilcox<-data.frame(random_wilcox)
random_wilcox1<-subset(random_wilcox,P.value<=1)
pathway1high_name<-subset(rownames(random_wilcox1),random_wilcox1[,4] == "High")
pathway2high_name<-subset(rownames(random_wilcox1),random_wilcox1[,4] == "Low")

dat_cluster<-cbind(TCGA_GSVA[,Cluster11],TCGA_GSVA[,Cluster22])
dat_cluster_pathway<-rbind(dat_cluster[pathway1high_name,],dat_cluster[pathway2high_name,])

common<-intersect(colnames(dat_cluster_pathway),rownames(pdata))
pdata<-pdata[common,]
pdata<-pdata[order(pdata[,4]),]
heatmap<-dat_cluster_pathway[,rownames(pdata)]

col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#00AFFF", "white", "#FC0E00")
)

library("RColorBrewer")
a<-names(table(pdata[,3]))
my_col<-brewer.pal(n = length(a), name = "Set3")
names(my_col)<-a
col_an<-columnAnnotation(Subtype = factor(pdata[,4]),
                         Age = factor(pdata[,1]),
                         Gender = factor(pdata[,2]),
                         Stage = factor(pdata[,3]), 
                         col = list(Subtype = c("High" = "#da191c", "Low" = "#2E9FDF"),
                                    Age = c("Old" = "lightcyan","Young" = "sienna1"),
                                    Gender=c("FEMALE" = "yellow","MALE" = "lightpink"),
                                    Stage=my_col),
                         annotation_legend_param = list(
                           title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                           labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                           grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")),
                         border = F,
                         annotation_name_side = "right",
                         annotation_name_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold")
)

ha11=rowAnnotation(Pvalue=row_anno_text(as.matrix(random_wilcox1[rownames(dat_cluster_pathway),5]), rot =0,location=unit(0, "mm")))

ha_cn1=Heatmap(heatmap,
        col = col_fun,
        show_column_names = F,
        top_annotation = col_an,
        row_title = "",
        column_title = "Hallmark",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        row_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        cluster_columns = F,
        cluster_rows = F,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8,fontfamily = "serif"),
        heatmap_legend_param = list(title = "NES",
                                    title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                                    labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                                    legend_height = unit(8, "cm"),
                                    grid_width = unit(0.8, "cm")),
        border = TRUE
)
ha_cn1A<-ha_cn1+ha11
draw(ha_cn1A)

pdf("TCGA_hallmark.pdf",,width = 10,height = 8)
draw(ha_cn1A)
dev.off()

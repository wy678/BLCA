rm(list=ls())

library("cowplot")
library(clusterProfiler)
library(enrichplot)
library(tidyverse)

setwd("xx\\xxxx\\xxxx\\xxxx\\xxxx")

load("xx/xxxx/xxxx/xxxx/TCGA_BLCA_edgeR.rda")

hall<-read.gmt("h.all.v7.5.1.symbols.gmt")
hall<-separate(data = hall, col = term, into = c("Hallmark", "term"), sep = "HALLMARK_")
hall<-hall[,c(2:3)]

ge = nrDEG_edgeR$logFC
names(ge) = rownames(nrDEG_edgeR)
ge = sort(ge,decreasing = T)
head(ge)

em1 <- GSEA(ge, TERM2GENE = hall,
            nPerm=1000,
            verbose=FALSE,by="fgsea",
            pAdjustMethod="BH",
            pvalueCutoff=1)

save(em1,file = "GSEA_Hallmark.rda")

a<-order(abs(em1$NES),decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(em1, geneSetID = a[i], title = em1$Description[a[i]])
}

pdf("GSEA_BLCA_Hallmark.pdf",width = 20,height = 10)
plot_grid(gseaplot[[1]],gseaplot[[2]],
          gseaplot[[3]],gseaplot[[4]],
          gseaplot[[5]],gseaplot[[6]],
          gseaplot[[7]],gseaplot[[8]],nrow = 2)
dev.off()

################################

library(ggpubr)
library(ggrepel)

res<-em1@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+
  scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+
  scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05

p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+labs(x="NES",y="-log10(P-value)",title = "Hallmark")+
  theme(plot.title = element_text(size=15,family="serif",color="black",face= "bold",hjust = 0.5),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
        legend.title=element_text(size=14,family="serif",color="black",face="bold"),
        legend.text=element_text(size=14,family="serif",color="black",face="bold")
  )
p2

ggsave("Hallmark_scatter.pdf",width=8,height=8,units="in")


##################################3

cell_gene<-read.csv("immune.csv",header = T)
cell_gene <- cell_gene[,c(6,1)]

em2 <- GSEA(ge, TERM2GENE = cell_gene,
            nPerm=1000,
            verbose=FALSE,by="fgsea",
            pAdjustMethod="BH",
            pvalueCutoff=1)

save(em2,file = "GSEA_Immport.rda")

a<-order(abs(em2$NES),decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(em2, geneSetID = a[i], title = em2$Description[a[i]])
}

pdf("GSEA_BLCA_ImmPort.pdf",width = 20,height = 10)
plot_grid(gseaplot[[1]],gseaplot[[2]],
          gseaplot[[3]],gseaplot[[4]],
          gseaplot[[5]],gseaplot[[6]],
          gseaplot[[7]],gseaplot[[8]],nrow = 2)
dev.off()

#############################

res<-em2@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+
  scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+
  scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05

p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+labs(x="NES",y="-log10(P-value)",title = "Immport")+
  theme(plot.title = element_text(size=15,family="serif",color="black",face= "bold",hjust = 0.5),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
        legend.title=element_text(size=14,family="serif",color="black",face="bold"),
        legend.text=element_text(size=14,family="serif",color="black",face="bold")
  )
p2

ggsave("Immport_scatter.pdf",width=8,height=8,units="in")

###########################################

metabolic<-read.table("metabolic_gene-114.txt",sep = "\t",fill = T)
metabolic<-metabolic[,c(2,1)]

em3 <- GSEA(ge, TERM2GENE = metabolic,
            nPerm=1000,
            verbose=FALSE,by="fgsea",
            pAdjustMethod="BH",
            pvalueCutoff=1)

save(em3,file = "GSEA_metabolic.rda")

a<-order(abs(em3$NES),decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(em3, geneSetID = a[i], title = em3$Description[a[i]])
}

pdf("GSEA_BLCA_metabolic.pdf",width = 20,height = 10)
plot_grid(gseaplot[[1]],gseaplot[[2]],
          gseaplot[[3]],gseaplot[[4]],
          gseaplot[[5]],gseaplot[[6]],
          gseaplot[[7]],gseaplot[[8]],nrow = 2)
dev.off()

################################3



#############################

res<-em3@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+
  scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+
  scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05

p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+labs(x="NES",y="-log10(P-value)",title = "Metabolic")+
  theme(plot.title = element_text(size=15,family="serif",color="black",face= "bold",hjust = 0.5),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
        legend.title=element_text(size=14,family="serif",color="black",face="bold"),
        legend.text=element_text(size=14,family="serif",color="black",face="bold")
  )
p2

ggsave("metabolic_scatter.pdf",width=8,height=8,units="in")


###########################################


immunemarker<-read.table("immunemarker.txt",sep = "\t",fill = T,header = T)
immunemarker<-immunemarker[,c(2,1)]

em4 <- GSEA(ge, TERM2GENE = immunemarker,
            nPerm=1000,
            verbose=FALSE,by="fgsea",
            pAdjustMethod="BH",
            pvalueCutoff=1)

save(em4,file = "GSEA_immunemarker.rda")

a<-order(abs(em4$NES),decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(em4, geneSetID = a[i], title = em4$Description[a[i]])
}

pdf("GSEA_BLCA_immunemarker.pdf",width = 20,height = 10)
plot_grid(gseaplot[[1]],gseaplot[[2]],
          gseaplot[[3]],gseaplot[[4]],
          gseaplot[[5]],gseaplot[[6]],
          gseaplot[[7]],gseaplot[[8]],nrow = 2)
dev.off()

################################3



#############################

res<-em4@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+
  scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+
  scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05

p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+labs(x="NES",y="-log10(P-value)",title = "Immune marker")+
  theme(plot.title = element_text(size=15,family="serif",color="black",face= "bold",hjust = 0.5),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
        legend.title=element_text(size=14,family="serif",color="black",face="bold"),
        legend.text=element_text(size=14,family="serif",color="black",face="bold")
  )
p2

ggsave("immunemarker_scatter.pdf",width=8,height=8,units="in")


###########################################
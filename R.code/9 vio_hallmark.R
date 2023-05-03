
rm(list=ls())

library(GEOquery)
library(ggplot2)
library(circlize)
library(reshape2)
library(ggpubr)

setwd("xx\\xxxx\\xxxx\\xxxx")

load("xx/xxxx/xxxx/xxxx/TCGA_GSVA.rda")
load("xx/xxxx/xxxx/xxxx/risk_data.rda")

TCGA_GSVA<-data.frame(t(TCGA_GSVA))
risk_data<-risk_data[rownames(TCGA_GSVA),]
risk_data$group<-Hmisc::capitalize(risk_data$group)
TCGA_GSVA$Subtype<-risk_data$group

plotA=TCGA_GSVA[,c(51,1:50)]
plotAA=melt(plotA,id="Subtype")
plotAA<-na.omit(plotAA)
colnames(plotAA)=c("Subtype","Immune","NES");
plotAA<-data.frame(plotAA)
###############################################
dataA<-plotAA

#col<-c("cadetblue1", "palevioletred1")
col<-c("cadetblue1", "orange")

p1<-ggplot(dataA,aes(x=Subtype,y=NES,fill=Subtype))+
  geom_violin()+
  geom_boxplot(width=0.2,notch = T)+
  facet_wrap(~Immune,scales = "free",ncol=5)+
  scale_fill_manual(values=c("lightgreen","lightgoldenrod2"))+guides(fill=guide_legend(title="Subtype"))

p2<-p1+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=11,family="serif",color="black",face="bold"),#y???̶????ݵ???
             axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
             legend.title=element_text(size=11,family="serif",color="black",face="bold"),
             legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
#p<-p+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3<-p2+labs(x="Subtype",y="NES")+theme(strip.text.x=element_text(size=10,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.x=1.2,label.y=3.1)

p5B<-p4+theme(panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),
              legend.position = "none")
p5B


#p1<-p1+theme(legend.title = element_text())
TCGA_GSVA_vio<-p5B
save(TCGA_GSVA_vio,file="TCGA_GSVA_vio.rda")
ggsave("TCGA_GSVA_vioplot.pdf",width=15,height=30,units="in")

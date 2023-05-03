
############################
rm(list=ls())

library(ggpubr)
library(ggplot2)
library(reshape2)

setwd("xx\\xxxx\\xxxx\\xxxxx\\xxx")

load("xx/xxxx/xxxx/xxxx/risk_data.rda")

load("xx/xxxx/xxxx/xxxxx/VAEN_GDSC.A.pred_TCGA.rda")
load("xx/xxxx/xxxx/xxxxx/xxxx/blue_MMGS.rda")

colnames(risk_data)[4]<-"Risk score"
risk_data[,4]<-stringr::str_to_title(risk_data[,4])
rownames(risk_data)<-substr(rownames(risk_data),1,15)
drug<-GDSC[which(GDSC$Sample %in% rownames(risk_data)),]
rownames(drug)<-drug$Sample

hub=rownames(module_gene_MMGS[which(module_gene_MMGS[,1]>0.5&module_gene_MMGS[,2]>0.4),])
hubdata<-drug[,hub]
a<-risk_data[rownames(drug),]
hubdata$Sample<-rownames(hubdata)
risk_data<-risk_data[rownames(hubdata),]
risk_data$Sample<-hubdata$Sample

datA<-melt(hubdata,id.vars = "Sample")
colnames(datA)<-c("Sample","Drug","Value")
mydata<-merge(risk_data,datA,id.vars = "Sample")
mydata<-na.omit(mydata)

col<-c("cadetblue1", "palevioletred1")
col<-c("cadetblue1", "orange")

p<-ggboxplot(mydata, x="Drug", y="Value", color = "Risk score",size=0.8,notch = T,
             palette =c("#00AFBB", "#E7B800"), add = "jitter")
p1<-p+stat_compare_means(label = "p.format",
                         aes(`Risk score`=`Risk score`),method="wilcox.test")
p2<-p1+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),#y????
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             legend.title=element_text(size=12,family="serif",color="black",face="bold"),
             legend.text=element_text(size=12,family="serif",color="black",face="bold"),
)
p3<-p2+theme(panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
             panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
             panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p4<-p3+labs(x=NULL,y=NULL)+facet_wrap(.~Drug,scales= "free",nrow = 2)+
  theme(strip.text.x=element_text(size=10,colour="black",
                                  face="bold",family="serif"))

p4

drug_BLCA<-p4
save(drug_BLCA,file="drug_BLCA.rda")
ggsave("drug_BLCA.pdf",width=8,height=8,units="in")

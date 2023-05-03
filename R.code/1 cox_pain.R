
########################################
rm(list =ls())
library(Biobase)
library(clusterProfiler)
library(HGNChelper)
library(reshape2)
library(tibble)
setwd("xx\\xxxx\\xxxx\\xxxx\\")
load("TCGA_BLCA.rda")
setwd("xx\\xxxx\\xxxx\\xxxx\\xxx\\")


exp<-exprs(TCGA_BLCA)
pain<-read.gmt("GOBP_SENSORY_PERCEPTION_OF_PAIN.v2022.1.Hs.gmt")
a<-checkGeneSymbols(pain[,2])
b<-checkGeneSymbols(rownames(exp))

library(limma)
c<-strsplit2(b[,3],split = " /// ")
b1<-cbind(b,c)[,c(1,4,5,6)]
b2<-melt(b1,id.var = "x")
b2<-na.omit(b2)
b2<-b2[-c(which(b2[,3] == "")),c(1,3)]

exp<-data.frame(exp[b2[,1],])
exp<-rownames_to_column(exp,var = "gene")
exp$gene<-b2[,2]

pain_exp<-exp[which(exp$gene %in% pain[,2]),]
pain_exp<-pain_exp[-which(rownames(pain_exp) ==  "13973"),]

rownames(pain_exp)<-pain_exp$gene
pain_exp_BLCA<-pain_exp[,-which(colnames(pain_exp) == "gene")]

save(pain_exp_BLCA,file = "pain_exp_BLCA.rda")



###############################################
rm(list =ls())
library(Biobase)
library(survminer)
library(survcomp)
setwd("xx\\xxxxxx\\xxxx\\xxxx\\")
load("TCGA_BLCA.rda")
setwd("xx\\xxxxxx\\xxxx\\xxxx\\xxx\\")
load("pain_exp_BLCA.rda")

pdata<-pData(TCGA_BLCA)
pcoxdata<-pdata[,c("OS.time","OS")]
colnames(pcoxdata)<-c("times","status")
colnames(pain_exp_BLCA)<-gsub("\\.","-",colnames(pain_exp_BLCA))
pain_exp_BLCA<-t(pain_exp_BLCA)

coxdata<-cbind(pcoxdata,pain_exp_BLCA)

cox_pain<-matrix(nrow=(ncol(coxdata)),ncol=6)
for (i in 3:(ncol(coxdata))) {
  Bcox<-coxph(Surv(times, status)~coxdata[,i],data=coxdata)
  summcph<-summary(Bcox)
  cox_pain[i,1]<-summcph$conf.int[1]
  cox_pain[i,2]<-summcph$conf.int[3]
  cox_pain[i,3]<-summcph$conf.int[4]
  cox_pain[i,4]<-as.matrix(summcph$logtest)[3]
  cox_pain[i,5]<-as.matrix(summcph$sctest)[3]
  cox_pain[i,6]<-summcph$coefficients[5]
  print(i)
}
rownames(cox_pain)<-colnames(coxdata)
cox_pain<-na.omit(cox_pain)
colnames(cox_pain)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
save(cox_pain,file = "cox_pain.rda")

#################


load("cox_pain.rda")

cox_pain<-data.frame(cox_pain[which(cox_pain[,6]<0.05),])

hr=c(sprintf("%0.2f",as.numeric(cox_pain$HR)))
hr1=format(hr,digits=3)
low=c(round(as.numeric(cox_pain$Lower.95),2))
low1=format(low,digits = 3)
low2=gsub(" ","", low1)
upp=c(round(as.numeric(cox_pain$Upper.95),2))
upp1=format(upp,digits = 3)
upp2=gsub(" ","", upp1)

pvalue=c(as.numeric(cox_pain$p_value))
value=format(pvalue,scientific = T,digits = 3)
HR=paste0(hr1,"(",low2,"-",upp2,")")
set<-rownames(cox_pain)
dat=cbind(c("Feature",set),c("HR (95% CI)",HR),c("P-value",value))
datA<-dat[-1,]
datA<-data.frame(datA)

HR_cox<-data.frame(cox_pain)
HR_cox$Group<-row.names(HR_cox)
HR_cox$var<-row.names(HR_cox)
HR_cox$CI<-datA$X2
HR_cox$P.value<-datA$X3
HR_cox<-HR_cox[order(HR_cox$HR),]
HR_cox$gene<-factor(rownames(HR_cox),levels = rownames(HR_cox))


####################################
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cowplot)

HR_cox1<-HR_cox

cbbPalette<-colorRampPalette(brewer.pal(12,"Set3"))(nrow(HR_cox1))

p1<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$gene,family="serif",color=c("black"),
            hjust=0,fontface="bold",inherit.aes = FALSE,size=4) +
  ggtitle("Pain gene")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.065))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()
  )

p1<-p1+scale_fill_manual(values=cbbPalette)
p1


p2<-ggplot(HR_cox1,aes(gene,HR,fill =gene)) +
  theme(panel.background = element_rect(fill='transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black',size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face="bold",colour = "black",hjust = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("") +
  #labs(title="95% confidence intervals")+
  ggtitle("95% confidence intervals")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))


p2<-p2+geom_errorbar(aes(ymin=Lower.95,ymax=Upper.95,color=gene), 
                     width=0.5,size = 0.8) +
  geom_point(aes(color=gene),shape=22,size=6)+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  geom_hline(aes(yintercept=median(HR_cox1$HR)),linetype='dashed')+
  theme(#y???̶????ݵ???
    axis.text.x=element_text(size=13,family="serif"),
  )

p2 

p3<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$P.value,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("P-value")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p3


p4<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$CI,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("HR (95% CI)")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p4

p<-p1+p2+p4+p3+plot_layout(widths=c(2,5,2.5,1.5))

forestplot_pain<-p
forestplot_pain

save(forestplot_pain,file = "forestplot_pain.rda")

pdf('forestplot_pain.pdf',width = 10,height = 6)
forestplot_pain
dev.off()

#######################

HR_cox2<-HR_cox[-which(rownames(HR_cox) == "OPRD1"),]

cbbPalette<-colorRampPalette(brewer.pal(12,"Set3"))(nrow(HR_cox2))

p1<-ggplot(HR_cox2,aes(HR_cox2$gene,HR_cox2$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox2$gene,family="serif",color=c("black"),
            hjust=0,fontface="bold",inherit.aes = FALSE,size=4) +
  ggtitle("Pain gene")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.065))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()
  )

p1<-p1+scale_fill_manual(values=cbbPalette)
p1


p2<-ggplot(HR_cox2,aes(gene,HR,fill =gene)) +
  theme(panel.background = element_rect(fill='transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black',size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face="bold",colour = "black",hjust = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("") +
  #labs(title="95% confidence intervals")+
  ggtitle("95% confidence intervals")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))


p2<-p2+geom_errorbar(aes(ymin=Lower.95,ymax=Upper.95,color=gene), 
                     width=0.5,size = 0.8) +
  geom_point(aes(color=gene),shape=22,size=6)+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  geom_hline(aes(yintercept=median(HR_cox2$HR)),linetype='dashed')+
  theme(#y???̶????ݵ???
    axis.text.x=element_text(size=13,family="serif"),
  )

p2 

p3<-ggplot(HR_cox2,aes(HR_cox2$gene,HR_cox2$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox2$P.value,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("P-value")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p3


p4<-ggplot(HR_cox2,aes(HR_cox2$gene,HR_cox2$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox2$CI,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("HR (95% CI)")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p4

p<-p1+p2+p4+p3+plot_layout(widths=c(2,5,2.5,1.5))

forestplot_pain16<-p
forestplot_pain16

save(forestplot_pain16,file = "forestplot_pain16.rda")

pdf('forestplot_pain16.pdf',width = 10,height = 6)
forestplot_pain16
dev.off()
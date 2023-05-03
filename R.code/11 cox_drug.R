
rm(list = ls())

library(survminer)
library(survcomp)
setwd("xx\\xxxx\\xxxx\\xxxxx\\xxx")

load("xx/xxxx/xxxx/xxxxx/VAEN_GDSC.A.pred_TCGA.rda")
load("xx/xxxx/xxxx/xxxxx/xxxx/blue_MMGS.rda")
load("xx/xxxx/xxxx/xxxx/risk_data.rda")

rownames(risk_data)<-substr(rownames(risk_data),1,15)
drug<-GDSC[which(GDSC$Sample %in% rownames(risk_data)),]
rownames(drug)<-drug$Sample

hub=rownames(module_gene_MMGS[which(module_gene_MMGS[,1]>0.5&module_gene_MMGS[,2]>0.4),])
hubdata<-drug[,hub]
a<-risk_data[rownames(drug),]
coxdata<-cbind(a,hubdata)

cox_drug<-matrix(nrow=(ncol(coxdata)),ncol=6)
for (i in 5:(ncol(coxdata))) {
  Bcox<-coxph(Surv(times, status)~coxdata[,i],data=coxdata)
  summcph<-summary(Bcox)
  cox_drug[i,1]<-summcph$conf.int[1]
  cox_drug[i,2]<-summcph$conf.int[3]
  cox_drug[i,3]<-summcph$conf.int[4]
  cox_drug[i,4]<-as.matrix(summcph$logtest)[3]
  cox_drug[i,5]<-as.matrix(summcph$sctest)[3]
  cox_drug[i,6]<-summcph$coefficients[5]
  print(i)
}
rownames(cox_drug)<-colnames(coxdata)
cox_drug<-na.omit(cox_drug)
colnames(cox_drug)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
save(cox_drug,file = "cox_drug.rda")

#################

load("cox_drug.rda")

cox_drug<-data.frame(cox_drug[which(cox_drug[,6]<0.05),])

hr=c(sprintf("%0.2f",as.numeric(cox_drug$HR)))
hr1=format(hr,digits=3)
low=c(round(as.numeric(cox_drug$Lower.95),2))
low1=format(low,digits = 3)
low2=gsub(" ","", low1)
upp=c(round(as.numeric(cox_drug$Upper.95),2))
upp1=format(upp,digits = 3)
upp2=gsub(" ","", upp1)

pvalue=c(as.numeric(cox_drug$p_value))
value=format(pvalue,scientific = T,digits = 3)
HR=paste0(hr1,"(",low2,"-",upp2,")")
set<-rownames(cox_drug)
dat=cbind(c("Feature",set),c("HR (95% CI)",HR),c("P-value",value))
datA<-dat[-1,]
datA<-data.frame(datA)

HR_cox<-data.frame(cox_drug)
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
  ggtitle("Drug")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.065))+
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

forestplot_drug<-p
forestplot_drug

save(forestplot_drug,file = "forestplot_drug.rda")

pdf('forestplot_drug.pdf',width = 10,height = 6)
forestplot_drug
dev.off()


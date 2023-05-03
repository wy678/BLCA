
rm(list=ls())

setwd("xx\\xxxx\\xxxx\\xxxxx\\xxxx")

library(reshape2)
library(tibble)
library(clusterProfiler)
library(HGNChelper)
library(Biobase)
library(survminer)
library(survcomp)
library(cowplot)

load("xx/xxx/xxxx/xxxxx/IMvigor210_clinical.rda")

exp<-exprs(IMvigor210)

pain<-read.gmt("xx\\xxxx\\xxxx\\xxxx\\GOBP_SENSORY_PERCEPTION_OF_PAIN.v2022.1.Hs.gmt")
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
pain_exp<-pain_exp[-which(rownames(pain_exp) ==  "17422"),]

rownames(pain_exp)<-pain_exp$gene
pain_exp_IM210<-pain_exp[,-which(colnames(pain_exp) == "gene")]

save(pain_exp_IM210,file = "pain_exp_IM210.rda")


pain_exp_IM210<-t(pain_exp_IM210)
pdata<-pData(IMvigor210)
pdata<-pdata[,c("os","censOS")]

coxdata<-cbind(pdata,pain_exp_IM210)

cox_pain<-matrix(nrow=(ncol(coxdata)),ncol=6)
for (i in 3:(ncol(coxdata))) {
  Bcox<-coxph(Surv(os, censOS)~coxdata[,i],data=coxdata)
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
cox_pain_IM210<-na.omit(cox_pain)
colnames(cox_pain_IM210)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
save(cox_pain_IM210,file = "cox_pain_IM210.rda")

####################################

cox_pain_IM210<-data.frame(cox_pain_IM210[which(cox_pain_IM210[,6]<0.05),])
pain_exp_IM210<-t(pain_exp_IM210)
pain_exp_IM2101<-as.matrix(pain_exp_IM210[rownames(cox_pain_IM210),])

data<-matrix(NA,nrow(pain_exp_IM2101),ncol(pain_exp_IM2101))
colnames(data)<-colnames(pain_exp_IM2101)
rownames(data)<-rownames(pain_exp_IM2101)

for (k in 1:nrow(pain_exp_IM2101)) {
  data[k,]<-pain_exp_IM2101[k,]*log(cox_pain_IM210[k,1])
}

data1<-data.frame(t(data))
data1$score<-apply(data1,1,sum)

risk_data<-pdata[rownames(data1),]
colnames(risk_data)<-c("times","status")
risk_data$score<-data1$score


sur.cut <- surv_cutpoint(risk_data,time = "times", event = "status",variables ="score",minprop=0.4)
sur.cat <- surv_categorize(sur.cut)
fit1 <- survfit(Surv(times, status)~score, data = sur.cat)
Bcox<-coxph(Surv(times, status)~score, data = sur.cat)
summcph<-summary(Bcox)
p1=ggsurvplot(fit1,legend.title="Subtype",
              main="Overall survival",
              xlab="Time (month)",
              ylab="Survival probability",
              title=paste0("Kaplan-Meier Curve for IMvigor210"),
              legend.labs=c("High","Low"),
              pval=F,conf.int = TRUE, 
              risk.table =T,
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(10, "bold", "black"),
                                      font.caption = c(10, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(12,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")),
              palette = c("#E7B800", "#2E9FDF") )
p2<-p1$plot+theme(plot.title = element_text(hjust = 0.5))
p3<-p2+annotate("text", x=5, y=0.25, label=paste0("P=",signif(summcph$sctest[3], 3)),color="black",family="serif",cex=5)
p4<-plot_grid(p3, p1$table, ncol = 1, align = 'v',rel_heights = c(3,1))
p4

pdf("IM210_pain_sur.pdf",width = 8,height = 6)
p4
dev.off()

risk_data$group<-sur.cat$score
risk_data_IM210<-risk_data
save(risk_data_IM210,file = "risk_data_IM210.rda")

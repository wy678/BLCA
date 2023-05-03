
###################################
rm(list=ls())
library("survival")
library("survminer")
library("cowplot")
setwd("x\\xxxxx\\xxx\\xxxx\\")
load("TCGA_BLCA.rda")
setwd("x\\xxxx\\xxxx\\xxxx\\xxxx\\")
load("cox_pain.rda")
load("pain_exp_BLCA.rda")

cox_pain1<-data.frame(cox_pain[which(cox_pain[,6]<0.05),])

pain_exp_BLCA1<-as.matrix(pain_exp_BLCA[rownames(cox_pain1),])

data<-matrix(NA,nrow(pain_exp_BLCA1),ncol(pain_exp_BLCA1))
colnames(data)<-colnames(pain_exp_BLCA1)
rownames(data)<-rownames(pain_exp_BLCA1)

for (k in 1:nrow(pain_exp_BLCA1)) {
  data[k,]<-pain_exp_BLCA1[k,]*log(cox_pain1[k,1])
}

data1<-data.frame(t(data))
data1$score<-apply(data1,1,sum)




pdata<-pData(TCGA_BLCA)
risk_data<-pdata[,c("OS.time","OS")]
colnames(risk_data)<-c("times","status")

rownames(data1)<-gsub("\\.","-",rownames(data1))
data1<-data1[rownames(risk_data),]

risk_data$score<-data1$score
risk_data$times<-risk_data$times/30

sur.cut <- surv_cutpoint(risk_data,time = "times", event = "status",variables ="score",minprop=0.4)
sur.cat <- surv_categorize(sur.cut)
fit1 <- survfit(Surv(times, status)~score, data = sur.cat)
Bcox<-coxph(Surv(times, status)~score, data = sur.cat)
summcph<-summary(Bcox)
p1=ggsurvplot(fit1,legend.title="Subtype",
              main="Overall survival",
              xlab="Time (month)",
              ylab="Survival probability",
              title=paste0("Kaplan-Meier Curve for TCGA cohort"),
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
p3<-p2+annotate("text", x=10, y=0.25, label=paste0("P=",signif(summcph$sctest[3], 3)),color="black",family="serif",cex=5)
p4<-plot_grid(p3, p1$table, ncol = 1, align = 'v',rel_heights = c(3,1))
p4

pdf("BLCA_pain_sur.pdf",width = 8,height = 6)
p4
dev.off()

risk_data$group<-sur.cat$score
save(risk_data,file = "risk_data.rda")


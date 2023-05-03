
#########################################
rm(list=ls())
setwd("x:\\xxxx\\xxx\\xxxxxx")

library(magrittr)
library(ggpubr)
library(Hmisc)
library(reshape2)

load("x:/xxx/xxx/xxx/BLCA_estimateScore.rda")
load("x:/xxx/xxx/xxx/risk_data.rda")

BLCA<-estimateScore[["TCGA-BLCA_FPKM"]]
BLCA<-BLCA[rownames(risk_data),]
BLCA$Risk<-capitalize(risk_data$group)
BLCA_plot<-melt(BLCA,id.vars = "Risk")
BLCA_plot<-na.omit(BLCA_plot)

BLCA_plot$value<-as.numeric(BLCA_plot$value)
CYT<-BLCA_plot[which(BLCA_plot$variable == "CYT"),]

outliers <-boxplot(CYT$value, plot=FALSE)$out
CYT<- CYT[-which(CYT$value %in% outliers),]

p1<-CYT %>% ggplot(aes(Risk,value,fill=Risk)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_violin()+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4,notch = T)+
  facet_grid(~variable,space="free",scale="free_x")+
  scale_fill_manual(values=c("#FC4E07","#00AFBB"))+
  stat_compare_means(aes(x = Risk,y = value),
                     method="wilcox.test",label = "p.format",label.x = 1.4)+
  scale_size_continuous(range=c(1,3))+
  labs(x=NULL,y=NULL)+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
        axis.text.y = element_text(color="black",size=10,face="bold"),
        axis.text.x = element_text(color="black",size=10,vjust=0.5,hjust = 1,
                                   face="bold",angle = 90),
        axis.line.x.top  = element_line(color="black"), 
        axis.text.x.top = element_blank(),
        axis.ticks.y.right=element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.x.top=element_blank(),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        text = element_text(size = 15,face =  "bold",family="serif"),
        strip.text.x = element_text(size=18,family="serif",color="black",face="bold"),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        legend.title=element_text(size=13,family="serif",color="black",face="bold"),
        legend.text=element_text(size=10,family="serif",color="black",face="bold")
  )+guides(x.sec="axis",y.sec = "axis")
p1

ggsave("TCGA_CYT.pdf",width = 6,height = 6)

others<-BLCA_plot[which(BLCA_plot$variable != "CYT"),]

p2<-others %>% ggplot(aes(Risk,value,fill=Risk)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_violin()+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4,notch = T)+
  facet_wrap(~variable,scale="free_y")+
  scale_fill_manual(values=c("#FC4E07","#00AFBB"))+
  stat_compare_means(aes(x = Risk,y = value),
                     method="wilcox.test",label = "p.format",label.x = 1.4)+
  scale_size_continuous(range=c(1,3))+
  labs(x=NULL,y=NULL)+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
        axis.text.y = element_text(color="black",size=10,face="bold"),
        axis.text.x = element_text(color="black",size=10,vjust=0.5,hjust = 1,
                                   face="bold",angle = 90),
        axis.line.x.top  = element_line(color="black"), 
        axis.text.x.top = element_blank(),
        axis.ticks.y.right=element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.x.top=element_blank(),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        text = element_text(size = 15,face =  "bold",family="serif"),
        strip.text.x = element_text(size=18,family="serif",color="black",face="bold"),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        legend.title=element_text(size=13,family="serif",color="black",face="bold"),
        legend.text=element_text(size=10,family="serif",color="black",face="bold")
  )+guides(x.sec="axis",y.sec = "axis")
p2
ggsave("TCGA_estimate.pdf",width = 10,height = 10)

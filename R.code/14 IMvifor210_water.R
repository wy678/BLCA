

rm(list = ls())

library(ggplot2)
library(ggpubr)
library(Hmisc)

setwd("xx\\xxxx\\xxxx\\xxxxx\\xxxx")

load("xx/xxxx/xxxx/xxxx/xxxx/risk_data_IM210.rda")
load("xx/xxxx/xxxx/xxxx/IMvigor210_clinical.rda")

pdata<-pData(IMvigor210)
score<-risk_data_IM210[rownames(pdata),]
score$binary<-pdata$binaryResponse

data<-na.omit(score)
data$group<-capitalize(data$group)
data$sample<-rownames(data)

box_data<-data[,c(3,5)]

p<-ggplot(data = box_data,aes(x = binary,y =score,fill = binary))+
  geom_boxplot()+
  scale_x_discrete(name =NULL) +
  scale_y_continuous(name = "Risk score")+
  scale_fill_manual(values = c("slateblue1","darkorange"))+
  theme_bw()+
  theme(legend.text = element_text(size=14),
        axis.title = element_text(size=15),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10),
        text = element_text(size = 12,face =  "bold",family="serif"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"))+ 
  stat_compare_means(aes(x = binary,y =score),
                     method="wilcox.test",label = "p.format",label.x = 1.4)
p

ggsave("IMvigor_boxplot.pdf",width = 6,height = 6)

#############################

mydata<-data[,c(4,5)]
colnames(mydata)<-c("Risk score","type")

chisqdata<-tidyr::unite(mydata, "res_risk",`Risk score`,type)
chisq<-chisq.test(matrix(table(chisqdata[,1]),nrow = 2,ncol = 2))

ggplot(mydata, aes(x = `Risk score`, fill = type)) +
  geom_bar(width = 0.8, position = "fill") + # ?ٷֱ???״ͼ
  annotate("text", label = paste0("Chisq.test p-value = ",signif(chisq$p.value,3)), x = 1.5, 
           y = 1.05, size = 5,family="serif")+
  scale_fill_manual(values=c("#da191c", "#2E9FDF")) +
  guides(fill=guide_legend(title = NULL)) +
  labs(x=NULL,y=NULL)+
  theme(legend.text = element_text(size = 15,face =  "bold",family="serif"),
        legend.title = element_text(size = 15,face =  "bold",family="serif"),
        axis.text.y =element_text(size = 10,face =  "bold",family="serif"),
        legend.position = "top",
        axis.text.x = element_text(size = 10,face =  "bold",family="serif"),
        axis.ticks=element_blank(),
        panel.background = element_blank())

ggsave("IMvigor_barplot.pdf",width = 6,height = 6)

############################

waterdata<-risk_data_IM210
waterdata$binary<-pdata$Best.Confirmed.Overall.Response

data<-waterdata[,c(3,5)]

data1<-data[which(data[,2]=="CR"|data[,2]=="PR"|data[,2]=="PD"|data[,2]=="SD"),]
data1$sample<-rownames(data1)
data1<-data1[order(data1$score),]
data1$sample<-factor(data1$sample,levels=as.vector(data1$sample))

library(wesanderson)

p<-ggplot(data=data1,aes(x=sample,y=score,fill=binary))+
  geom_bar(position=position_dodge(), stat="identity")+
  scale_fill_manual(values = wes_palette('Darjeeling1'))+
  labs(title="IMvigor210",
       x="Anti-PD-L1 Response", y = "Risk score", fill = "")+
  theme_bw()+
  theme(legend.position='right',
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="serif",face="bold",size=12,color="Black"),
        legend.text = element_text(family="serif",colour="black", size = 12, face = "bold"),
        axis.title = element_text(family="serif",size = 15, face="bold", color = "Black"),
        plot.title = element_text(family="serif",size=18,hjust = 0.5,lineheight=.8, face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p


ggsave(p,filename = "IMvigor210-waterplot.pdf",width = 10,height = 6)


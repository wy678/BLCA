

############################
rm(list=ls())

library(ggpubr)
library(ggplot2)
library(reshape2)

setwd("xx\\xxxx\\xxxx\\xxxx")

load("xx/xxxx/xxxx/xxxx/risk_data.rda")
a<-read.table("Mat_BLCA.tsv",sep = "\t",header = T)
colnames(risk_data)[4]<-"Risk score"
risk_data[,4]<-stringr::str_to_title(risk_data[,4])

bar<-a[,c(1,3)]
colnames(bar)<-c("Sample","type")
risk_data$Sample<-rownames(risk_data)
risk_data<-risk_data[,-c(1,2,3)]
mydata<-merge(risk_data,bar,id.vars = "Sample")

mydata$type<-apply(data.frame(mydata[,3]),2,function(x){
  ifelse(x=="False","Non-responder","Responser")})
mydata<-na.omit(mydata)


chisqdata<-tidyr::unite(mydata, "res_risk",`Risk score`,type)
chisq<-chisq.test(matrix(table(chisqdata[,2]),nrow = 2,ncol = 2))

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

ggsave("TIDE_bar.pdf",width = 6,height = 8)

############################################

box<-a[,c(1,4,10:14)]
datA<-melt(box,by="X")
colnames(datA)<-c("Sample","Gene","Value")
mydata<-merge(risk_data,datA,id.vars = "Sample")
mydata<-na.omit(mydata)

col<-c("cadetblue1", "palevioletred1")
col<-c("cadetblue1", "orange")

p<-ggboxplot(mydata, x="Gene", y="Value", color = "Risk score",size=0.8,notch = T,
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
p4<-p3+labs(x=NULL,y=NULL)+facet_wrap(.~Gene,scales= "free",nrow = 2)+
  theme(strip.text.x=element_text(size=10,colour="black",
                                  face="bold",family="serif"))

p4

TIDE_BLCA<-p4
save(TIDE_BLCA,file="TIDE_BLCA.rda")
ggsave("TIDE_BLCA.pdf",width=8,height=8,units="in")

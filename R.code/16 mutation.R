#################################################################
rm(list=ls())

setwd("xx\\xxxx\\xxxx\\xxxxx")

library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(maftools)
library(survminer)


load("xx/xxxx/xxx/xxxx/maftools-BLCA.rda")
load("xx/xxxx/xxx/xxxx/risk_data.rda")

#########################################


dat=laml@data
sample1=unique(as.character(dat$Tumor_Sample_Barcode))

a<-laml@variants.per.sample
a$x1<-substr(a$Tumor_Sample_Barcode,1,16)
rownames(risk_data)<-substr(rownames(risk_data),1,16)
a<-a[which(a$x1 %in% rownames(risk_data)),]
risk_data<-risk_data[a$x1,]
rownames(risk_data)<-a$Tumor_Sample_Barcode

colnames(risk_data)[4]<-"Subtype"
clu_high<-rownames(risk_data[which(risk_data$Subtype=="high"),])
clu_low<-rownames(risk_data[which(risk_data$Subtype=="low"),])

clu_high_maf=subsetMaf(laml,tsb = clu_high)
clu_low_maf=subsetMaf(laml,tsb = clu_low)

save(clu_high_maf,clu_low_maf,file="maf_Risk.rda")

dat1=laml@data
gene=unique(as.character(dat1$Hugo_Symbol))
spname=unique(as.character(dat1$Tumor_Sample_Barcode))
matt=matrix(nrow = length(gene),ncol = length(spname))
colnames(matt)=spname
rownames(matt)=gene
times=c()
for (i in 1:length(gene)) {
  a=dat1[which(as.character(dat1$Hugo_Symbol)==gene[i]),]
  times[i]=nrow(a)
  for (j in 1:nrow(a)) {
    matt[i,as.character(a$Tumor_Sample_Barcode)[j]]=as.character(a$Variant_Classification)[j]
  }
}
save(matt,file="BLCAmatt.rda")


matt_clu_high=matt[,clu_high]
matt_clu_low=matt[,clu_low]
times_clu_high=c()
for (i in 1:nrow(matt_clu_high)) {
  a=matt_clu_high[i,-which(is.na(matt_clu_high[i,]))]
  times_clu_high[i]=length(a)
}

times_clu_low=c()
for (i in 1:nrow(matt_clu_low)) {
  b=matt_clu_low[i,-which(is.na(matt_clu_low[i,]))]
  times_clu_low[i]=length(b)
}

gene_clu_high=rownames(matt_clu_high)
names(times_clu_high)=gene_clu_high
gene_clu_high_top=names(sort(times_clu_high,decreasing = T)[1:20])
matt_clu_high1=matt_clu_high[gene_clu_high_top,]
matt_clu_high1[which(is.na(matt_clu_high1))]=""

gene_clu_low=rownames(matt_clu_low)
names(times_clu_low)=gene_clu_low
gene_clu_low_top=names(sort(times_clu_low,decreasing = T)[1:20])
matt_clu_low1=matt_clu_low[gene_clu_low_top,]
matt_clu_low1[which(is.na(matt_clu_low1))]=""

mutates=unique(as.character(dat1$Variant_Classification))
type1=c()
type2=c()
for (i in 1:length(mutates)) {
  if(length(which(matt_clu_high1==mutates[i]))>0)
    type1=c(type1,mutates[i])
  if(length(which(matt_clu_low1==mutates[i]))>0)
    type2=c(type2,mutates[i])
}

type=union(type1,type2)
#type=c("Missense_Mutation","Frame_Shift_Del","Splice_Site",
#       "In_Frame_Del","Nonsense_Mutation","Frame_Shift_Ins",
#       "In_Frame_Ins","Nonstop_Mutation")

type=c("Frame_Shift_Ins","Missense_Mutation","Splice_Site",
       "Nonsense_Mutation","Frame_Shift_Del","In_Frame_Del",
       "In_Frame_Ins", "Translation_Start_Site","Nonstop_Mutation")

#col=brewer.pal(12,"Set3")
#col=brewer.pal(length(type),"Set3")
col= c("Frame_Shift_Ins"="#E6AB02",
       "Missense_Mutation" = "#1B9E77",
       "Splice_Site" = "#7570B3",
       "Nonsense_Mutation"="#66A61E",
       "Frame_Shift_Del" = "#D95F02",
       "In_Frame_Del"="#E7298A", 
       "In_Frame_Ins"="#A6761D",
       "Translation_Start_Site"="#CCEBC5",
       "Nonstop_Mutation"="#FFED6F" 
)

names(col)=type

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "white", col = NA))
  },
  "Missense_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"),gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  "Frame_Shift_Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h*0.5, gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  "Splice_Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"),h-unit(0.5, "mm"), gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  "In_Frame_Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  "Nonsense_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  "Frame_Shift_Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  "In_Frame_Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  
  "Nonstop_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  },
  "Translation_Start_Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.35, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  }
)


colnames(risk_data)[4]<-"Risk score"
risk_data[,4]<-stringr::str_to_title(risk_data[,4])
#????ע????Ϣ
ha_clu_high<-HeatmapAnnotation('Risk score'=risk_data[clu_high,"Risk score"],show_annotation_name = F,show_legend = F,
                               col = list('Risk score' = c("High" =  "red")),
                               annotation_legend_param =list('Risk score'=list(title="Risk score",
                                                                               title_gp = gpar(fontsize = 11,fontface = "bold",fontfamily="serif"),
                                                                               labels_gp = gpar(fontsize = 8,fontface = "bold",fontfamily="serif")))
                               
)

ha_clu_low<-HeatmapAnnotation('Risk score'=risk_data[clu_low,"Risk score"],show_annotation_name = F,show_legend = F,
                              col = list('Risk score' = c("Low" =  "blue")),
                              annotation_legend_param =list('Risk score'=list(title="Risk score",
                                                                              title_gp = gpar(fontsize = 11,fontface = "bold",fontfamily="serif"),
                                                                              labels_gp = gpar(fontsize = 8,fontface = "bold",fontfamily="serif")))
                              
)



#ָ?????????͵ı?ǩ?????????е????Ͷ?Ӧ
heatmap_legend_param <- list(title = "Alternations",at =type,labels = type,
                             title_gp = gpar(fontsize = 11,fontface = "bold",fontfamily="serif"),
                             labels_gp = gpar(fontsize = 8,fontface = "bold",fontfamily="serif")
)

clu_high_p=oncoPrint(matt_clu_high1,alter_fun = alter_fun, col =col,row_names_side = "left", #????????
                     show_heatmap_legend = F,
                     remove_empty_columns=TRUE,
                     remove_empty_rows=TRUE,
                     pct_side = "right",pct_gp = gpar(fontsize = 10,fontface = "bold",fontfamily="serif"),
                     bottom_annotation = ha_clu_high, #ע????Ϣ?ڵײ?
                     heatmap_legend_param = heatmap_legend_param
                     #,right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(bar_width = 0, height = NULL,axis = F))
)

clu_high_p


clu_low_p=oncoPrint(matt_clu_low1,alter_fun = alter_fun, col =col,row_names_side = "left", #????????
                    show_heatmap_legend=F,
                    remove_empty_columns=TRUE,
                    remove_empty_rows=TRUE,
                    pct_side="right",pct_gp = gpar(fontsize = 10,fontface = "bold",fontfamily="serif"),
                    bottom_annotation=ha_clu_low, #ע????Ϣ?ڵײ?
                    heatmap_legend_param=heatmap_legend_param
                    #,right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(bar_width = 0,border = F,axis = F))
                    
)

clu_low_p

pp=clu_high_p+clu_low_p
lgd_list = list(
  Legend(at=c("clu_high","clu_low"),labels = c("High","Low"),labels_gp = gpar(fontsize=9,fontface="bold",fontfamily="serif"),
         title="Risk score",title_gp=gpar(fontsize = 10,fontface="bold",fontfamily="serif"),
         legend_gp = gpar(fill = c("red","blue"))),
  Legend(title = "Alternations",title_gp = gpar(fontsize=10,fontface = "bold",fontfamily="serif"),
         at =type,labels = type,
         labels_gp = gpar(fontsize = 9,fontface = "bold",fontfamily="serif"),
         legend_gp = gpar(fill =col))
  #legend_gp = gpar(fill = brewer.pal(9,"Set3")))
)

draw(pp, ht_gap = unit(1, "mm"), annotation_legend_list=lgd_list,heatmap_legend_side="bottom",
     column_title="OncoPrint for BLCA cohort",column_title_gp=gpar(fontsize = 18,fontface="bold",fontfamily="serif"))

pdf("BLCA_mutation_OncoPrint_Risk.pdf",height=6.5,width=14.5)
draw(pp, ht_gap = unit(2.5, "mm"), annotation_legend_list = lgd_list,
     column_title="OncoPrint for BLCA cohort",column_title_gp=gpar(fontsize = 18,fontface = "bold",fontfamily="serif"))
dev.off()







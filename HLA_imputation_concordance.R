##Check the imputation concordance between CookHLA and HIBAG
##Concordance was calculated by summing the number of matching alleles between CookHLA and HIBAG across participants 
##and dividing by twice the number of participants
##https://www.nature.com/articles/s42003-023-05496-5#Fig6 (concordance definition )

setwd("/Users/jingjing/Documents/R_files/HLA_Concordance")
library(rio)
library(tidyverse)
library(ggplot2)
##Compare allele frequency----
PSI_HIBAG<-import("PSI_HIBAG_Prefit_imputation.frq",format = "tsv")
PSI_T1DGC<-import("PSI_T1DGC_REF.MHC.HLA_IMPUTATION_OUT.bmarker.FRQ.frq",format="tsv")

##subset classical HLA alleles 
HLA_PSI_HIBAG <- PSI_HIBAG[grepl("HLA",PSI_HIBAG$ID),]
HLA_PSI_T1DGC<- PSI_T1DGC[grepl("HLA",PSI_T1DGC$SNP),]
intersect(HLA_PSI_HIBAG$ID,HLA_PSI_T1DGC$SNP)

##Keep the shared HLA alleles
HLA_Com<-HLA_PSI_HIBAG %>% left_join(HLA_PSI_T1DGC[,2:6],by=c("ID"="SNP")) %>% na.omit()
HLA_Com$group<- str_sub(HLA_Com$ID, end=-7)
##make a scatter plot to see the consistency

pdf(file="Allele_frequency.pdf")
F<-ggplot(data=HLA_Com,aes(x=HLA_Com$ALT_FREQS,y=HLA_Com$MAF,color = group))+
  geom_point()+
  geom_smooth(alpha=0.1,method='lm',formula= y~x,aes(fill=group))+
  theme_prism(palette = "warm_pastels",
             base_size = 11, 
             base_fontface = "plain", 
             border = TRUE) + 
  scale_colour_prism(palette = "warm_pastels") + 
  xlab("Allele_frequency_HIBAG")+
  ylab("Allele_frequency_CookHLA")+
  stat_cor()


F
dev.off()
##Separate plot for each HLA classical allele----
library(ggpubr)
pdf(file="Separated_alleles.pdf")
P<-ggscatter(HLA_Com, 
             size=1.5,
             x = "ALT_FREQS", 
             y ="MAF" ,
             color="group",
             add = "reg.line",
             alpha=0.4) +
  facet_wrap(~group) +
  stat_cor(label.y = 0.4,size=3,digits = 3,aes(color=group),show.legend=FALSE)+
  geom_smooth(alpha=0.2,method='lm',formula= y~x,aes(fill=group,color=group,alpha=0.2),lwd=0.2)+
  theme_prism(
              palette = "warm_pastels",
              base_size = 11, 
              base_fontface = "plain", 
              border = TRUE) + 
  scale_colour_prism(palette = "warm_pastels") + 
  xlab("Allele_frequency_HIBAG")+
  ylab("Allele_frequency_CookHLA")
P
dev.off()
##Calculate correlation of allele frequency between two methods----
correlate <- HLA_Com %>%
  group_by(group) %>% 
 summarise(r=cor(ALT_FREQS, MAF))


##Sum number of matched alleles----
PSI_HIBAG_allele<- import("PSI_HIBAG_Prefit_imputation.imgt3560.2field.chped.subset.chped",format = "tsv")
PSI_T1DGC_allele<- import("PSI_T1DGC_REF.MHC.HLA_IMPUTATION_OUT.imgt3560.2field.chped",format = "tsv")

##HLA_A matching----
PSI_HLA_com<- PSI_HIBAG_allele[,c(1,7,8)]
PSI_HLA_com_1<-PSI_T1DGC_allele[,c(1,7,8)]
PSI_HLA_A_matching<- PSI_HLA_com_1 %>% left_join(PSI_HLA_com,by=c("FID"="FID"))
colnames(PSI_HLA_A_matching)<-c("FID","CookHLA_A1","CookHLA_A2","HIBAG_A1","HIBAG_A2")

HLA_A_unmatched<-PSI_HLA_A_matching[PSI_HLA_A_matching$CookHLA_A1!=PSI_HLA_A_matching$HIBAG_A1&PSI_HLA_A_matching$CookHLA_A1!=PSI_HLA_A_matching$HIBAG_A2,]
HLA_A_unmatched_2<-PSI_HLA_A_matching[PSI_HLA_A_matching$CookHLA_A2!=PSI_HLA_A_matching$HIBAG_A2&PSI_HLA_A_matching$CookHLA_A2!=PSI_HLA_A_matching$HIBAG_A1,]
HLA_A_unmatched_total<-rbind(HLA_A_unmatched,HLA_A_unmatched_2)
write.table(HLA_A_unmatched_total,file = "HLA_A_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")
##HLA_A unmatched 63/(3465*2)=0.00909
##concordance rate= 1-63/(3465*2)=99.09% for HLA_A

##HLA_B matching----
PSI_HLA_B_HIBAG<- PSI_HIBAG_allele[,c(1,9,10)]
PSI_HLA_B_T1DGC<-PSI_T1DGC_allele[,c(1,9,10)]
PSI_HLA_B_matching<- PSI_HLA_B_T1DGC %>% left_join(PSI_HLA_B_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_B_matching)<-c("FID","CookHLA_B1","CookHLA_B2","HIBAG_B1","HIBAG_B2")

HLA_B_unmatched <- PSI_HLA_B_matching[PSI_HLA_B_matching$CookHLA_B1!=PSI_HLA_B_matching$HIBAG_B1&PSI_HLA_B_matching$CookHLA_B1!=PSI_HLA_B_matching$HIBAG_B2,]
HLA_B_unmatched_2<- PSI_HLA_B_matching[PSI_HLA_B_matching$CookHLA_B2!=PSI_HLA_B_matching$HIBAG_B2&PSI_HLA_B_matching$CookHLA_B2!=PSI_HLA_B_matching$HIBAG_B1,]
HLA_B_unmatched_total<-rbind(HLA_B_unmatched,HLA_B_unmatched_2)
HLA_B_unmatched_total<-distinct(HLA_B_unmatched_total)
write.table(HLA_B_unmatched_total,file = "HLA_B_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")
##HLA_B unmatched 114/(3465*2) = 0.0164
##concordance rate= 1-114/(3465*2) = 98.35%

##HLA_C matching----
PSI_HLA_C_HIBAG<-PSI_HIBAG_allele[,c(1,11,12)]
PSI_HLA_C_T1DGC<-PSI_T1DGC_allele[,c(1,11,12)]
PSI_HLA_C_matching<- PSI_HLA_C_T1DGC %>% left_join(PSI_HLA_C_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_C_matching)<-c("FID","CookHLA_C1","CookHLA_C2","HIBAG_C1","HIBAG_C2")

HLA_C_unmatched <- PSI_HLA_C_matching[PSI_HLA_C_matching$CookHLA_C1!=PSI_HLA_C_matching$HIBAG_C1&PSI_HLA_C_matching$CookHLA_C1!=PSI_HLA_C_matching$HIBAG_C2,]
HLA_C_unmatched_2<- PSI_HLA_C_matching[PSI_HLA_C_matching$CookHLA_C2!=PSI_HLA_C_matching$HIBAG_C2&PSI_HLA_C_matching$CookHLA_C2!=PSI_HLA_C_matching$HIBAG_C1,]
HLA_C_unmatched_total<-rbind(HLA_C_unmatched,HLA_C_unmatched_2)
HLA_C_unmatched_total<-distinct(HLA_C_unmatched_total)
write.table(HLA_C_unmatched_total,file = "HLA_C_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")
##HLA_C unmatched 31/(3465*2) =0.00447
##HLA_C concordance rate= 1-31/(3465*2)=99.55%

##HLA_DPB1 matching----
PSI_HLA_DPB1_HIBAG<- PSI_HIBAG_allele[,c(1,13,14)]
PSI_HLA_DPB1_T1DGC<- PSI_T1DGC_allele[,c(1,15,16)]
PSI_HLA_DPB1_matching<- PSI_HLA_DPB1_T1DGC %>% left_join(PSI_HLA_DPB1_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_DPB1_matching)<-c("FID","CookHLA_DPB1_1","CookHLA_DPB1_2","HIBAG_DPB1_1","HIBAG_DPB1_2")

HLA_DPB1_unmatched<-PSI_HLA_DPB1_matching[PSI_HLA_DPB1_matching$CookHLA_DPB1_1!=PSI_HLA_DPB1_matching$HIBAG_DPB1_1 &
                                          PSI_HLA_DPB1_matching$CookHLA_DPB1_1!=PSI_HLA_DPB1_matching$HIBAG_DPB1_2,]

HLA_DPB1_unmatched_2<-PSI_HLA_DPB1_matching[PSI_HLA_DPB1_matching$CookHLA_DPB1_2!=PSI_HLA_DPB1_matching$HIBAG_DPB1_2&
                                            PSI_HLA_DPB1_matching$CookHLA_DPB1_2!=PSI_HLA_DPB1_matching$HIBAG_DPB1_1,]


HLA_DPB1_unmatched_total<-rbind(HLA_DPB1_unmatched,HLA_DPB1_unmatched_2) 
write.table(HLA_DPB1_unmatched_total,file = "HLA_DPB1_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")
##HLA_DPB1_unmatched: 315/(3465*2)= 0.04545455
##HLA_DPB1_concordance rate: 1-315/(3465*2)= 95.45%

#Find unique individuals
HLA_DPB1_unmatched_total_distinct<-distinct(HLA_DPB1_unmatched_total)

##HLA_DQA1 matching----
PSI_HLA_DQA1_HIBAG<-PSI_HIBAG_allele[,c(1,15,16)]
PSI_HLA_DQA1_T1DGC<- PSI_T1DGC_allele[,c(1,17,18)]
PSI_HLA_DQA1_matching<- PSI_HLA_DQA1_T1DGC %>% left_join(PSI_HLA_DQA1_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_DQA1_matching)<-c("FID","CookHLA_DQA1_1","CookHLA_DQA1_2","HIBAG_DQA1_1","HIBAG_DQA1_2")

HLA_DQA1_unmatched<-PSI_HLA_DQA1_matching[PSI_HLA_DQA1_matching$CookHLA_DQA1_1!=PSI_HLA_DQA1_matching$HIBAG_DQA1_1 &
                                          PSI_HLA_DQA1_matching$CookHLA_DQA1_1!=PSI_HLA_DQA1_matching$HIBAG_DQA1_2,]

HLA_DQA1_unmatched_2<-PSI_HLA_DQA1_matching[PSI_HLA_DQA1_matching$CookHLA_DQA1_2!=PSI_HLA_DQA1_matching$HIBAG_DQA1_2&
                                            PSI_HLA_DQA1_matching$CookHLA_DQA1_2!=PSI_HLA_DQA1_matching$HIBAG_DQA1_1,]

HLA_DQA1_unmatched_total<-rbind(HLA_DQA1_unmatched,HLA_DQA1_unmatched_2)
write.table(HLA_DQA1_unmatched_total,file = "HLA_DQA1_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")
##HLA_DQA1_unmatched:(462+852)/(3465*2)=0.1896
##HLA_DQA1_concordance rate:1-(462+852)/(3465*2)=81.04%

##HLA_DQB1 matching----
PSI_HLA_DQB1_HIBAG<-PSI_HIBAG_allele[,c(1,17,18)]
PSI_HLA_DQB1_T1DGC<- PSI_T1DGC_allele[,c(1,19,20)]
PSI_HLA_DQB1_matching<- PSI_HLA_DQB1_T1DGC %>% left_join(PSI_HLA_DQB1_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_DQB1_matching)<-c("FID","CookHLA_DQB1_1","CookHLA_DQB1_2","HIBAG_DQB1_1","HIBAG_DQB1_2")

HLA_DQB1_unmatched<-PSI_HLA_DQB1_matching[PSI_HLA_DQB1_matching$CookHLA_DQB1_1!=PSI_HLA_DQB1_matching$HIBAG_DQB1_1 &
                                            PSI_HLA_DQB1_matching$CookHLA_DQB1_1!=PSI_HLA_DQB1_matching$HIBAG_DQB1_2,]

HLA_DQB1_unmatched_2<-PSI_HLA_DQB1_matching[PSI_HLA_DQB1_matching$CookHLA_DQB1_2!=PSI_HLA_DQB1_matching$HIBAG_DQB1_2&
                                              PSI_HLA_DQB1_matching$CookHLA_DQB1_2!=PSI_HLA_DQB1_matching$HIBAG_DQB1_1,]

HLA_DQB1_unmatched_total<-rbind(HLA_DQB1_unmatched,HLA_DQB1_unmatched_2)
write.table(HLA_DQB1_unmatched_total,file = "HLA_DQB1_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")

##HLA_DQB1_unmatched: 29/(3465*2)=0.00418
##HLA_DQB1_concordance: 1-29/(3465*2)=99.58%

##HLA_DRB1 matching----
PSI_HLA_DRB1_HIBAG<-PSI_HIBAG_allele[,c(1,19,20)]
PSI_HLA_DRB1_T1DGC<-PSI_T1DGC_allele[,c(1,21,22)]
PSI_HLA_DRB1_matching<- PSI_HLA_DRB1_T1DGC %>% left_join(PSI_HLA_DRB1_HIBAG,by=c("FID"="FID"))
colnames(PSI_HLA_DRB1_matching)<-c("FID","CookHLA_DRB1_1","CookHLA_DRB1_2","HIBAG_DRB1_1","HIBAG_DRB1_2")

HLA_DRB1_unmatched<-PSI_HLA_DRB1_matching[PSI_HLA_DRB1_matching$CookHLA_DRB1_1!=PSI_HLA_DRB1_matching$HIBAG_DRB1_1 &
                                            PSI_HLA_DRB1_matching$CookHLA_DRB1_1!=PSI_HLA_DRB1_matching$HIBAG_DRB1_2,]

HLA_DRB1_unmatched_2<-PSI_HLA_DRB1_matching[PSI_HLA_DRB1_matching$CookHLA_DRB1_2!=PSI_HLA_DRB1_matching$HIBAG_DRB1_2&
                                              PSI_HLA_DRB1_matching$CookHLA_DRB1_2!=PSI_HLA_DRB1_matching$HIBAG_DRB1_1,]

HLA_DRB1_unmatched_total<-rbind(HLA_DRB1_unmatched,HLA_DRB1_unmatched_2)
write.table(HLA_DRB1_unmatched_total,file = "HLA_DRB1_unmatched_total.csv",quote = F, row.names = F, col.names = T,sep="\t")

##HLA_DRB1_unmatched:179/(3465*2) =0.0258
##HLA_DRB1_concordance:1-179/(3465*2)=97.42%


##make a barplot to show the concordance rate across all the shared classical alleles----
HLA_names<- c("HLA-A","HLA-B","HLA-C","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRB1")
Concordance<- c(0.9909,0.9835,0.9955,0.9545,0.8104,0.9958,0.9742)
df<-data.frame(HLA_names,Concordance,row.names=NULL)
library(ggprism)
pdf(file = "HLA_Concordance.pdf")
P<-ggplot(df,aes(x=HLA_names,y=Concordance))+
  geom_bar(stat = "identity",aes(fill=factor(HLA_names)))+ 
  theme_prism(palette = "warm_pastels",
              axis_text_angle = 45,
              base_size = 11, 
              base_fontface = "plain", 
              base_family = "Arial",
              border = TRUE) + 
  scale_colour_prism(palette = "warm_pastels") + 
  scale_fill_prism(palette = "warm_pastels")+
  xlab("") 

P
dev.off()


##Prepare files for the circos plot 
##keep the chromosome information
hg38<-data.frame(chromosome=c(1:22),
                 size=c(248956422,
                        242193529,
                        198295559,
                        190214555,
                        181538259,
                        170805979,
                        159345973,
                        145148636,
                        148394717,
                        143797422,
                        145086622,
                        143275309,
                        114364328,
                        107043718,
                        101991189,
                        90338345,
                        83257441,
                        80373285,
                        58617616,
                        64444167,
                        46709983,
                        50818468))

write.table(hg38,file="hg38.csv",sep="\t",quote= F, col.names=T,row.names=F)
##First layer ever/never
library(tidyverse)
library(readr)
library(readxl)
IBD_ever<-read_xlsx("Meta_IBD_ever_lead_Annotation.xlsx")
IBD_ever<-IBD_ever[c(2,4,5,11)]
colnames(IBD_ever)[1]<-"gSNP"
IBD_ever$`Gene context` <-str_replace_all(IBD_ever$`Gene context`,"\\]","")
IBD_ever$`Gene context` <-str_replace_all(IBD_ever$`Gene context`,"–––","/")
IBD_ever$`Gene context` <-str_replace_all(IBD_ever$`Gene context`,",","/")
IBD_ever$`Gene context` <-str_replace_all(IBD_ever$`Gene context`,"/RAB9BP1","RAB9BP1")
write.table(IBD_ever,file="cir_IBD_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_ever")
CD_ever<-read_xlsx("Meta_CD_ever_lead_Annotation.xlsx")
CD_ever<-CD_ever[c(2,4,5,11)]
CD_ever$`Gene context` <-str_replace_all(CD_ever$`Gene context`,"––\\[\\]\\––","/")
CD_ever$`Gene context` <-str_replace_all(CD_ever$`Gene context`,"–\\[\\]\\-","/")
CD_ever$`Gene context` <-str_replace_all(CD_ever$`Gene context`,",","/")
write.table(CD_ever,file="cir_CD_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_ever")
UC_ever<-read_xlsx("Meta_UC_ever_lead_Annotation.xlsx")
UC_ever<-UC_ever[c(2,4,5,11)]
UC_ever$`Gene context` <-str_replace_all(UC_ever$`Gene context`,"BCL10-AS1BCL10","BCL10-AS1/BCL10")
write.table(UC_ever,file="cir_UC_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

##current
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_current")
IBD_current<-read_xlsx("Meta_IBD_current_lead_Annotation.xlsx")
IBD_current<-IBD_current[c(2,4,5,11)]
IBD_current$`Gene context` <-str_replace_all(IBD_current$`Gene context`,",","/")
write.table(IBD_current,file="cir_IBD_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_current")
CD_current<-read_xlsx("Meta_CD_current_lead_Annotation.xlsx")
CD_current<-CD_current[c(2,4,5,11)]
CD_current$`Gene context` <-str_replace_all(CD_current$`Gene context`,"––\\[\\]\\––","/")
CD_current$`Gene context` <-str_replace_all(CD_current$`Gene context`,",","")
write.table(CD_current,file="cir_CD_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_current")
UC_current<-read_xlsx("Meta_UC_current_lead_Annotation.xlsx")
UC_current<-UC_current[c(2,4,5,11)]
UC_current$`Gene context` <-str_replace_all(UC_current$`Gene context`,"––\\[\\]\\––","/")
UC_current$`Gene context` <-str_replace_all(UC_current$`Gene context`,"-\\[\\]\\–","/")
UC_current$`Gene context` <-str_replace_all(UC_current$`Gene context`,",","")
write.table(UC_current,file="cir_UC_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

##former
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_former/")
IBD_former<-read_xlsx("Meta_IBD_former_lead_Annotation.xlsx")
IBD_former<-IBD_former[c(2,4,5,11)]
IBD_former$`Gene context` <-str_replace_all(IBD_former$`Gene context`,",","/")
write.table(IBD_former,file="cir_IBD_former.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_former/")
CD_former<-read_xlsx("Meta_CD_former_lead_Annotation.xlsx")
CD_former<-CD_former[c(2,4,5,11)]
CD_former$`Gene context` <-str_replace_all(CD_former$`Gene context`,"B3GAT1/","B3GAT1")
write.table(CD_former,file="cir_CD_former.csv",sep="\t",quote= F, col.names=T,row.names=F)

setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_former/")
UC_former<-read_xlsx("Meta_UC_former_lead_Annotation.xlsx")
UC_former<-UC_former[c(2,4,5,11)]
UC_former$`Gene context` <-str_replace_all(UC_former$`Gene context`,"\\]","")
write.table(UC_former,file="cir_UC_former.csv",sep="\t",quote= F, col.names=T,row.names=F)

##add yadav
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_ever/Yadav/")
published<-read_xlsx("Yadav_Anno.xlsx")
published<-published[c(2,4,5,11)]
published$Group<-"Known"
published$`Gene context` <-str_replace_all(published$`Gene context`,"\\([^()]{0,}\\)","")
published$`Gene context` <-str_replace_all(published$`Gene context`,"TSHR/","TSHR")

##combine all the gene names 
combine<- bind_rows(IBD_ever,
                    IBD_current,
                    IBD_former,
                    CD_ever,
                    CD_current,
                    CD_former,
                    UC_current,
                    UC_ever,
                    UC_former,
                    published)
combine[is.na(combine$Group),]$Group <-"New"
##subset the first chromosome
combine1<-combine[combine$Chr =="1",]
combine1$No<-c(1:52)
##remove the variants with same position or gene name 
combine1[duplicated(combine1$`Gene context`) ==TRUE,]
combine1[duplicated(combine1$Pos) ==TRUE,]
combine1<-combine1[-c(4,14,24,43,44,23),]
##second
combine2<-combine[combine$Chr =="2",]
combine2$No<-c(1:41)
combine2[duplicated(combine2$`Gene context`) ==TRUE,]
combine2[duplicated(combine2$Pos) ==TRUE,]
combine2<-combine2[-c(19,20,22,23),]
##third
combine3<-combine[combine$Chr =="3",]
combine3$No<-c(1:54)
combine3[duplicated(combine3$`Gene context`) ==TRUE,]
combine3[duplicated(combine3$Pos) ==TRUE,]
combine3<-combine3[-c(9,27,29,31,32,48),]
##fourth
combine4<-combine[combine$Chr =="4",]
combine4$No<-c(1:49)
combine4[duplicated(combine4$`Gene context`) ==TRUE,]
combine4[duplicated(combine4$Pos) ==TRUE,]
combine4<-combine4[-c(8,22,26,30,33,34,37,45,47,20),]
##fifth
combine5<-combine[combine$Chr =="5",]
combine5$No<-c(1:49)
combine5[duplicated(combine5$`Gene context`) ==TRUE,]
combine5[duplicated(combine5$Pos) ==TRUE,]
combine5<-combine5[-c(7,18,24,26,30,37,43),]
##sixth
combine6<-combine[combine$Chr =="6",]
combine6$No<-c(1:54)
combine6[duplicated(combine6$`Gene context`) ==TRUE,]
combine6[duplicated(combine6$Pos) ==TRUE,]
combine6<-combine6[-c(1,29,2,43,10,37,14,41,53),]
##seventh
combine7<-combine[combine$Chr =="7",]
combine7$No<-c(1:46)
combine7[duplicated(combine7$`Gene context`) ==TRUE,]
combine7[duplicated(combine7$Pos) ==TRUE,]
combine7<-combine7[-c(2,5,25,32,36,37,42,22),]
##eighth
combine8<-combine[combine$Chr =="8",]
combine8$No<-c(1:35)
combine8[duplicated(combine8$`Gene context`) ==TRUE,]
combine8[duplicated(combine8$Pos) ==TRUE,]
combine8<-combine8[-c(20,21,27,32,33),]
##nineth
combine9<-combine[combine$Chr =="9",]
combine9$No<-c(1:25)
combine9[duplicated(combine9$`Gene context`) ==TRUE,]
combine9[duplicated(combine9$Pos) ==TRUE,]
combine9<-combine9[-c(17,18),]
##tenth
combine10<-combine[combine$Chr =="10",]
combine10$No<-c(1:24)
combine10[duplicated(combine10$`Gene context`) ==TRUE,]
combine10[duplicated(combine10$Pos) ==TRUE,]
combine10<-combine10[-c(11,15,17,21),]
##eleventh
combine11<-combine[combine$Chr =="11",]
combine11$No<-c(1:34)
combine11[duplicated(combine11$`Gene context`) ==TRUE,]
combine11[duplicated(combine11$Pos) ==TRUE,]
combine11<-combine11[-c(22,26,28,30,33,14),]
##twelve
combine12<-combine[combine$Chr =="12",]
combine12$No<-c(1:29)
combine12[duplicated(combine12$`Gene context`) ==TRUE,]
combine12[duplicated(combine12$Pos) ==TRUE,]
combine12<-combine12[-c(14,16,26),]
##thirteenth
combine13<-combine[combine$Chr =="13",]
combine13$No<-c(1:17)
combine13[duplicated(combine13$`Gene context`) ==TRUE,]
combine13[duplicated(combine13$Pos) ==TRUE,]
combine13<-combine13[-c(8,5),]
##fourteenth
combine14<-combine[combine$Chr =="14",]
combine14$No<-c(1:16)
combine14[duplicated(combine14$`Gene context`) ==TRUE,]
combine14[duplicated(combine14$Pos) ==TRUE,]
combine14<-combine14[-c(2,13),]
##fifteenth
combine15<-combine[combine$Chr =="15",]
combine15$No<-c(1:20)
combine15[duplicated(combine15$`Gene context`) ==TRUE,]
combine15[duplicated(combine15$Pos) ==TRUE,]
combine15<-combine15[-c(12,18),]
##sisteenth
combine16<-combine[combine$Chr =="16",]
combine16$No<-c(1:23)
combine16[duplicated(combine16$`Gene context`) ==TRUE,]
combine16[duplicated(combine16$Pos) ==TRUE,]
combine16<-combine16[-c(11,15,16,18,20),]
##seventeenth
combine17<-combine[combine$Chr =="17",]
combine17$No<-c(1:21)
combine17[duplicated(combine17$`Gene context`) ==TRUE,]
combine17[duplicated(combine17$Pos) ==TRUE,]
combine17<-combine17[-c(19,20),]
##eighteenth
combine18<-combine[combine$Chr =="18",]
combine18$No<-c(1:17)
combine18[duplicated(combine18$`Gene context`) ==TRUE,]
combine18[duplicated(combine18$Pos) ==TRUE,]
combine18<-combine18[-c(10,11,12,13,14),]
##nineteenth
combine19<-combine[combine$Chr =="19",]
combine19$No<-c(1:9)
combine19[duplicated(combine19$`Gene context`) ==TRUE,]
combine19[duplicated(combine19$Pos) ==TRUE,]
combine19<-combine19[-c(7),]
#twenty
combine20<-combine[combine$Chr =="20",]
combine20$No<-c(1:9)
combine20[duplicated(combine20$`Gene context`) ==TRUE,]
combine20[duplicated(combine20$Pos) ==TRUE,]
combine20<-combine20[-c(7,9),]
##twenty-one
combine21<-combine[combine$Chr =="21",]
combine21$No<-c(1:7)
combine21[duplicated(combine21$`Gene context`) ==TRUE,]
combine21[duplicated(combine21$Pos) ==TRUE,]
combine21<-combine21[-c(4),]
##twenty two
combine22<-combine[combine$Chr =="22",]
combine22$No<-c(1:11)
combine22[duplicated(combine22$`Gene context`) ==TRUE,]
combine22[duplicated(combine22$Pos) ==TRUE,]
Cir_label<- bind_rows(combine1,combine2,
                      combine3,combine4,
                      combine5,combine6,
                      combine7,combine8,
                      combine9,combine10,
                      combine11,combine12,
                      combine13,combine14,
                      combine15,combine16,
                      combine17,combine18,
                      combine19,combine20,
                      combine21,combine22)
write.table(Cir_label,file="Cir_label.csv",sep="\t",quote= F, col.names=T,row.names=F)

published$Outcome<- c(rep("IBD_former",3),"IBD_ever",rep("IBD_current",3),"IBD_ever","IBD_current","IBD_current","IBD_current","IBD_former","IBD_current","IBD_former","IBD_ever_former","IBD_current","IBD_former","IBD_current","IBD_current","IBD_ever_current","IBD_ever_former","IBD_ever_current","IBD_current","IBD_former","IBD_former",rep("UC_current",9),"UC_ever","UC_former","UC_current","UC_current","UC_current","UC_current","UC_current","UC_current","UC_current","UC_current","UC_current","UC_current","CD_former","CD_former","CD_ever","CD_ever","CD_former","CD_former","CD_former","CD_former","CD_former","CD_former","CD_current","CD_current","CD_ever_current","CD_ever_current","CD_current","CD_ever_current","CD_former","UC_current")

Yadav_IBD_ever<-published[grep("IBD_ever",published$Outcome),]
write.table(Yadav_IBD_ever,file="Yadav_IBD_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_IBD_current<-published %>% filter(grepl("IBD",Outcome)) %>% filter(grepl("current",Outcome))
write.table(Yadav_IBD_current,file="Yadav_IBD_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_IBD_former<-published %>% filter(grepl("IBD",Outcome)) %>% filter(grepl("former",Outcome))
write.table(Yadav_IBD_former,file="Yadav_IBD_former.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_CD_ever<-published %>% filter(grepl("CD",Outcome)) %>% filter(grepl("ever",Outcome))
write.table(Yadav_CD_ever,file="Yadav_CD_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_CD_current<-published %>% filter(grepl("CD",Outcome)) %>% filter(grepl("current",Outcome))
write.table(Yadav_CD_current,file="Yadav_CD_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_CD_former<-published %>% filter(grepl("CD",Outcome)) %>% filter(grepl("former",Outcome))
write.table(Yadav_CD_former,file="Yadav_CD_former.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_UC_ever<-published %>% filter(grepl("UC",Outcome)) %>% filter(grepl("ever",Outcome))
write.table(Yadav_UC_ever,file="Yadav_UC_ever.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_UC_current<-published %>% filter(grepl("UC",Outcome)) %>% filter(grepl("current",Outcome))
write.table(Yadav_UC_current,file="Yadav_UC_current.csv",sep="\t",quote= F, col.names=T,row.names=F)

Yadav_UC_former<-published %>% filter(grepl("UC",Outcome)) %>% filter(grepl("former",Outcome))
write.table(Yadav_UC_former,file="Yadav_UC_former.csv",sep="\t",quote= F, col.names=T,row.names=F)




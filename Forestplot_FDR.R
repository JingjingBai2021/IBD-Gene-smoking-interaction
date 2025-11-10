library(SNPannotator)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(forestploter)
install.packages("psych")
library(psych)
detach(package:ggplot2, unload = T)
install.packages('inflection',dependencies=TRUE)
library("inflection")
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_ever/")

##read FDR 
IBD_ever_FDR <-readr::read_table("FDR.0.1")
IBD_ever_FDR
##identify the inflection point (FDR=0.288) 
d<-density(IBD_ever_FDR$FDR)
plot(d)
cc=check_curve(d$x,d$y)
ipese=ese(d$x,d$y,cc$index);ipese
ipbese=bese(d$x,d$y,cc$index);ipbese
ibede=bede(d$x,d$y,cc$index);ibede

pdf("IBD_ever_fdr_plot.pdf")
IBD_ever_FDR_plot
dev.off()
class(IBD_ever_FDR)
IBD_ever_inflection_plot<-IBD_ever_FDR %>%
  filter(IBD_ever_FDR$P<0.05) %>%
  ggplot(aes(x=P,y=FDR)) +
  geom_point(color="lightgrey", alpha=0.5) +
  ggtitle("IBD Ever/Never | FDR inflection=0.28") +
  geom_hline(yintercept=0.28, linetype="dashed", color = "#69b3a2")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
            
##read lead SNPs
Meta_IBD_ever_lead <- Meta_IBD_ever_lead[,c(1,4,3,5)]

##merge Lead snps with summary statistics
colnames(Meta_IBD_ever)
Meta_IBD_ever_lead <- Meta_IBD_ever_lead %>% left_join(Meta_IBD_ever[,c(3,4,7,9,11,16,17,19)],by =c("SNP"="SNP"))

##check missing value and consistency
which(is.na(Meta_IBD_ever_lead))
which(Meta_IBD_ever_lead$P.x!=Meta_IBD_ever_lead$P)

##adding FDR
Meta_IBD_ever_lead<- Meta_IBD_ever_lead %>% left_join(IBD_ever_FDR,by=c("SNP"="SNP"))

##remove the unnecessary ones
Meta_IBD_ever_lead$P.x<-NULL
Meta_IBD_ever_lead$P.y<-NULL

##order the columns 
Meta_IBD_ever_lead<-arrange(Meta_IBD_ever_lead,CHR,BP)

##combine the snp position columns 
Meta_IBD_ever_lead$"Chr:Pos" <- paste(Meta_IBD_ever_lead$CHR,Meta_IBD_ever_lead$BP,sep = ":")
Meta_IBD_ever_lead<-Meta_IBD_ever_lead[,-c(1,2)]
colnames(Meta_IBD_ever_lead)

##arrange the order of columns
Meta_IBD_ever_lead <- Meta_IBD_ever_lead[c(10,2,1,3,7,8,9,4,5,6)]

##annotate the lead_snps gene
db <- "1000GENOMES:phase_3:EUR"
server<- "https://rest.ensembl.org"
IBD_ever_lead_ann<- annotate(Meta_IBD_ever_lead$SNP,
                             server,db, 'Meta_IBD_ever_lead_Annotation.xlsx', 
                             LDlist = FALSE, 
                             cadd = FALSE, 
                             geneNames.file = '/Users/jingjing/Documents/R_files/UMCG/Data_base/Gene_Names_Ensembl_104_GRCh38.rds',
                             regulatoryType.file = 
                               '/Users/jingjing/Documents/R_files/UMCG/Data_base/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.rds')

IBD_ever_lead_ann<- readxl::read_xlsx("Meta_IBD_ever_lead_Annotation.xlsx")
Meta_IBD_ever_lead <- Meta_IBD_ever_lead %>% left_join(IBD_ever_lead_ann[,c(2,11)], by =c("SNP"="gSNP"))

##prepare forest plot
Meta_IBD_ever_lead$" " <- paste(rep(" ", 30), collapse = " ")
colnames(Meta_IBD_ever_lead)
Meta_IBD_ever_lead <- Meta_IBD_ever_lead[c(1,11,2,3,12,4,5,6,7,8)]
Meta_IBD_ever_lead$Q <-round(Meta_IBD_ever_lead$Q,digits = 3)
Meta_IBD_ever_lead$FDR <-round(Meta_IBD_ever_lead$FDR,digits = 3)
colnames(Meta_IBD_ever_lead)[3] <- "ALT"
Meta_IBD_ever_lead <- Meta_IBD_ever_lead %>% left_join(Meta_IBD_ever[,c(3,16,17)], by = c("SNP"="SNP"))
View(Meta_IBD_ever_lead)

png(filename="Meta_IBD_ever_lead.png", width = 18, height = 24,res = 600,units = "in")
library(extrafont)
extrafont::font_import()

tm<-forest_theme(base_family="Arial")
forest(Meta_IBD_ever_lead[,c(1,4,3,2,5,7,8,9,10)],
       est = Meta_IBD_ever_lead$OR,
       lower = Meta_IBD_ever_lead$Lower, 
       upper = Meta_IBD_ever_lead$Upper,
       ci_column = 5,
       ref_line = 1,
       arrow_lab = NULL,
       xlim = c(0, 8),
       ticks_at = c(0, 1, 2, 4, 6, 8))
       #footnote = "Brackets indicate SNP position.\nThe length of dashes indicate distance between SNP and flanking genes(-,>1Kb; –,>10 Kb; ––,>100 Kb)")

dev.off()
writexl::write_xlsx(Meta_IBD_ever_lead[,c(1,4,3,2,5,7,8,9,10)],"Meta_IBD_ever_lead.xlsx")
write.table(IBD_ever_lead, file ="Meta_IBD_ever_lead.csv",quote= F, col.names= T, row.names = F,sep ="\t")




setwd("/Users/jingjing/Documents/R_files/Hetero_CD_UC_swap/Ever_CD_UC/")
##Hetero_CD_UC_swap folder is the next version of opposite variants by swapping the major allele as effect allele
##extract the variants both are significant in CD an UC (P<0.01)----
library(readr)
library(readxl)
Meta_ever_CD<- readr::read_table("CD.meta")
Meta_ever_UC<-readr::read_table("UC.meta")
Overlapping_ever<- Reduce(intersect,
                          list(Meta_ever_CD[Meta_ever_CD$P<0.01,]$SNP,
                               Meta_ever_UC[Meta_ever_UC$P<0.01,]$SNP))
Overlapping_ever <- as.data.frame(Overlapping_ever)
colnames(Overlapping_ever)[1]<- "SNP"

##Add CD effect size----
library(tidyverse)
Overlapping_ever <- Overlapping_ever %>% 
                    left_join(Meta_ever_CD[,c(1:4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_ever)[5]<- "CD_P"
colnames(Overlapping_ever)[6]<- "CD_OR"
colnames(Overlapping_ever)[7]<- "CD_Q"

##Add UC effect size----
Overlapping_ever <- Overlapping_ever %>% left_join(Meta_ever_UC[,c(3,4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_ever)[9]<- "UC_P"
colnames(Overlapping_ever)[10]<- "UC_OR"

##Extract the overlapped variants those have opposite effect in CD and UC----
Opposite_ever<-Overlapping_ever[Overlapping_ever$CD_OR>1.0 & Overlapping_ever$UC_OR< 1.0,]
Opposite_ever2<-Overlapping_ever[Overlapping_ever$CD_OR<1.0 & Overlapping_ever$UC_OR> 1.0,]
Opposite_ever<-rbind(Opposite_ever,Opposite_ever2)

##Add position and ALT allele and 95% CI from summary statistics----
##CD
Meta_ever_CD$SE <- abs((log(Meta_ever_CD$OR))/qnorm(Meta_ever_CD$P/2))
Meta_ever_CD$"Lower" <- round(exp(log(Meta_ever_CD$OR)-1.96*Meta_ever_CD$SE),digits = 2)
Meta_ever_CD$"Upper" <- round(exp(log(Meta_ever_CD$OR)+1.96*Meta_ever_CD$SE),digits = 2)
Meta_ever_CD$"95% CI" <- paste(format(Meta_ever_CD$"Lower",digits=2),format(Meta_ever_CD$"Upper",digits=2), sep = "-")
Meta_ever_CD$"OR (95% CI)" <- paste0(format(Meta_ever_CD$OR,digits = 2)," (", Meta_ever_CD$`95% CI`, ")")
Meta_ever_CD$"OR (95% CI)"<- str_replace_all(Meta_ever_CD$"OR (95% CI)"," ","")
Meta_ever_CD$"95% CI" <-NULL

##UC
Meta_ever_UC$SE <- abs((log(Meta_ever_UC$OR))/qnorm(Meta_ever_UC$P/2))
Meta_ever_UC$"Lower" <- round(exp(log(Meta_ever_UC$OR)-1.96*Meta_ever_UC$SE),digits = 2)
Meta_ever_UC$"Upper" <- round(exp(log(Meta_ever_UC$OR)+1.96*Meta_ever_UC$SE),digits = 2)
Meta_ever_UC$"95% CI" <- paste(format(Meta_ever_UC$"Lower",digits=2),format(Meta_ever_UC$"Upper",digits=2), sep = "-")
Meta_ever_UC$"OR (95% CI)" <- paste0(format(Meta_ever_UC$OR,digits = 2)," (", Meta_ever_UC$`95% CI`, ")")
Meta_ever_UC$"OR (95% CI)"<- str_replace_all(Meta_ever_UC$"OR (95% CI)"," ","")
Meta_ever_UC$"95% CI" <-NULL

##add to 95%CI opposite variants
Opposite_ever<- Opposite_ever %>% left_join(Meta_ever_CD[,c(3,16,17,18)], by= c("SNP"="SNP"))
colnames(Opposite_ever)
Opposite_ever<-Opposite_ever[c(1,2,3,4,5,6,7,12,13,14,9,10,11)]
colnames(Opposite_ever)[8] <- "CD_Lower"
colnames(Opposite_ever)[9] <- "CD_Upper"
colnames(Opposite_ever)[10] <- "CD_OR(95%CI)"

Opposite_ever<- Opposite_ever %>% left_join(Meta_ever_UC[,c(3,16,17,18)], by= c("SNP"="SNP"))
colnames(Opposite_ever)[14] <- "UC_Lower"
colnames(Opposite_ever)[15] <- "UC_Upper"
colnames(Opposite_ever)[16] <- "UC_OR(95%CI)"
colnames(Opposite_ever)[4]<-"ALT"
##remove the missing value lines
Opposite_ever <- drop_na(Opposite_ever)

##save the files for obtaining independent variants
write.table(Opposite_ever[,c(1,2,3,4,10)],"Opposite_ever.csv",sep ="\t", quote=F, col.names=T,row.names=F)
write.table(Opposite_current[,c(1,2,3,4,10)],"Opposite_current.csv",sep ="\t", quote=F, col.names=T,row.names=F)
write.table(Opposite_former[,c(1,2,3,4,10)],"Opposite_former.csv",sep ="\t", quote=F, col.names=T,row.names=F)

##Make scatter plot within ever vs never contrast----
library(geneplotter)
##ever/never
##Select the ratio of odds ratio larger than 10 or smaller than 0.1
Opposite_ever$"Ratio of odds ratio"<-Opposite_ever$CD_OR/Opposite_ever$UC_OR
Highlight_ever<-subset(Opposite_ever,
                       Opposite_ever$"Ratio of odds ratio">10 |Opposite_ever$"Ratio of odds ratio"<0.1)
##annotate the highlighted snps
db <-"1000GENOMES:phase_3:EUR"
server<- "https://rest.ensembl.org"
Highlight_ever_density <- SNPannotator::annotate (Highlight_ever$SNP,
                                server,db, 'Highlight_ever_density_anno.xlsx', 
                                LDlist = FALSE, 
                                cadd = FALSE, 
                                geneNames.file = '/Users/jingjing/Documents/R_files/UMCG/Data_base/Gene_Names_Ensembl_104_GRCh38.rds',
                                regulatoryType.file ='/Users/jingjing/Documents/R_files/UMCG/Data_base/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.rds')

Highlight_ever_density_anno <-read_xlsx("Highlight_ever_density_anno.xlsx")
Highlight_ever <- Highlight_ever %>% left_join(Highlight_ever_density_anno[,c(2,11)], by =c("SNP"="gSNP"))
##keep the several snps mapped to the same gene
Highlight_ever<- subset(Highlight_ever,
                        Highlight_ever$`Ratio of odds ratio` <0.07|Highlight_ever$`Ratio of odds ratio`>10.382)
##make scatter plots
Lables<-c(0.0,1.0,2.0,3.0,4.0,5.0)
Lables<-as.numeric(Lables)
Lables<-sprintf(Lables,fmt='%#.1f')

pdf(file="Hetero_ever_CD_UC_new_2.pdf",width=6, height=7)
plot(Opposite_ever$CD_OR,Opposite_ever$UC_OR, 
     xlim=c(0,5),
     ylim=c(0,5),
     col=densCols(Opposite_ever$CD_OR,Opposite_ever$UC_OR),
     pch=20,      
     main = "Ever/Never smoking",
     xlab="CD",       
     ylab="UC",      
     axes = FALSE,
     cex.main = 1)
axis(2,at=c(0.0,1.0,2.0,3.0,4.0,5.0),tick=T,labels = Lables)
axis(1,at=c(0.0,1.0,2.0,3.0,4.0,5.0),tick=T,labels = Lables)
abline(coef = c(0,1),col="grey", lty=2)
text(Highlight_ever$CD_OR,
     Highlight_ever$UC_OR,
     Highlight_ever$`Gene context`,
     cex=0.65,     
     col="#164863",    
     offset =0.5,      
     font = 4)
points(Highlight_ever$CD_OR,
       Highlight_ever$UC_OR,
       col="#B51B75",    
       cex=0.8,    
       bg="#E1AFD1",    
       pch=21)      
dev.off()
##Making the density plot----
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(ggplot2)
pdf(file="Hetero_CD_UC_density_ever_new.pdf",width=6, height=7)
ggplot(data=Opposite_ever) +
  geom_density(aes(x=UC_OR), alpha=.4,col="black",fill="#F9E2AF") +
  geom_density(aes(x=CD_OR), alpha=.4,col="black",fill="#009FBD") +
  labs(title = "Ever/Never smoking", x= "OR", y= "Density")+
  theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color="black"))+
  scale_x_continuous(limits = c(0,5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) +
  annotate("text", x =1.8, y = 1.0, label = "CD",fontface="bold",col="#009FBD") +
  annotate("text", x = 0.4, y = 1.2, label = "UC",fontface="bold",col="#DE8F5F")

dev.off()

length(Opposite_ever[Opposite_ever$CD_OR >1.0,]$SNP)
length(Opposite_ever[Opposite_ever$UC_OR >1.0,]$SNP)


##Forest plot all contrast----
##Keep independent variants r2 <0.1 for forest plot----
##plink.clumped file was generated based on CD_statistics from UMCG cohort
Opposite_inde_ever<- readr::read_table(file = "/Users/jingjing/Downloads/plink_mac_20230116/ever/plink.clumped")
Opposite_inde_ever<-Opposite_inde_ever[,c(1,4,3,5)]
Opposite_inde_ever<- Opposite_inde_ever %>% left_join(Opposite_ever[,c(1,4:16)], by=c("SNP"="SNP"))
Opposite_inde_ever$P<-NULL

##Annotate the snps for forest plots----
##ever
library(readxl)
SNPs<-Opposite_inde_ever$SNP
Opposite_ever_anno <- annotate (SNPs,
                                server,db, 'Opposite_ever_anno.xlsx', 
                                LDlist = FALSE, 
                                cadd = FALSE, 
                                geneNames.file = '/Users/jingjing/Documents/R_files/UMCG/Data_base/Gene_Names_Ensembl_104_GRCh38.rds',
                                regulatoryType.file ='/Users/jingjing/Documents/R_files/UMCG/Data_base/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.rds')

Opposite_ever_anno<-read_xlsx("Opposite_ever_anno.xlsx")
Opposite_inde_ever <- Opposite_inde_ever %>% left_join(Opposite_ever_anno[,c(2,11)],by =c("SNP"="gSNP"))
Opposite_inde_ever<-Opposite_inde_ever[c(1:3,17,4:16)] 
Opposite_inde_ever<-arrange(Opposite_inde_ever,CHR,BP)
Opposite_inde_ever$"Chr:Pos" <- paste0(Opposite_inde_ever$CHR,":",Opposite_inde_ever$BP)
Opposite_inde_ever<-Opposite_inde_ever[c(24,3:23)] 
Opposite_inde_ever$P<- format(Opposite_inde_ever$P,scientific=T, digits=2)
Opposite_inde_ever$UC_P<- format(Opposite_inde_ever$UC_P,scientific=T, digits=2)

Opposite_inde_ever$CD<-NULL
colnames(Opposite_inde_ever)[20]<-""
colnames(Opposite_inde_ever)[15]<-"Pᵤ"
##make forest plot ever/never
library(forestploter)
Opposite_inde_ever$`CD` <- paste(rep(" ", 30), collapse = " ")
Opposite_inde_ever$"UC" <- paste(rep(" ", 30), collapse = " ")
tm<- forest_theme(base_size = 10,
                  refline_lty = "solid",
                  ci_pch = c(15,18),
                  ci_col = c("#009FBD", "#D21312"),
                  legend_name = "Group",
                  legend_value = c("CD", "UC"),
                  legend_position = "right",
                  vertline_lty = c("dotted"),
                  vertline_col = c("#bababa"))

tiff('Opposite_ever_forest_ALL.png', res = 600, width = 15, height = 30, units = "in")
unload(package = "ggplot2")
forest(Opposite_inde_ever[,c(1, 2, 3, 20, 10, 15)],
       est = list(Opposite_inde_ever$CD_OR,
                  Opposite_inde_ever$UC_OR),
       lower = list(Opposite_inde_ever$CD_Lower,
                    Opposite_inde_ever$UC_Lower), 
       upper = list(Opposite_inde_ever$CD_Upper,
                    Opposite_inde_ever$UC_Upper),
       ci_column = c(4),
       ref_line = 1,
       vert_line = c(0.5,1.5),
       xlim = c(0,12),
       nudge_y = 0.2,
       theme = tm)
dev.off()

Opposite_inde_ever_1<-Opposite_inde_ever[c(1:50),]
Opposite_inde_ever_2<-Opposite_inde_ever[c(51:101),]
png('Opposite_ever_forest_1.png', res = 600, width = 15, height = 17, units = "in")

forest(Opposite_inde_ever_1[,c(1, 2, 3, 20, 10, 15)],
       est = list(Opposite_inde_ever_1$CD_OR,
                  Opposite_inde_ever_1$UC_OR),
       lower = list(Opposite_inde_ever_1$CD_Lower,
                    Opposite_inde_ever_1$UC_Lower), 
       upper = list(Opposite_inde_ever_1$CD_Upper,
                    Opposite_inde_ever_1$UC_Upper),
       ci_column = c(4),
       ref_line = 1,
       vert_line = c(0.5,1.5),
       xlim = c(0,12),
       nudge_y = 0.2,
       theme = tm)
dev.off()





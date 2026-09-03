
##We added biological_use a new variate----
IBD_inflam_interaction<- read_delim("/Users/jingjing/Documents/R_files/RNA-seq/IBD_Inflam.cis_qtl_top_assoc.txt")
IBD_inflam_interaction_p_adj_0.05<- IBD_inflam_interaction %>% filter(IBD_inflam_interaction$pval_adj_bh<0.05) ##genome wide significant signals
View(IBD_inflam_interaction_p_adj_0.05) ## After correcting for biological use, there is not a major differences in the GxS eqtl in the results. 


filtered_data_qtl_bed_IBD_inflammation<-read_delim(file = "filtered_data_qtl_bed_IBD_inflammation.bed")
View(filtered_data_qtl_bed_IBD_inflammation)

##Check imputed results, Variants with maf>5% missing rate<5%.
IBD_inflam_interaction_imputed_v2<- read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/IBD_Inflam.cis_qtl_top_assoc.txt")
IBD_inflam_interaction_p_adj_0.05_imputed_v2<- IBD_inflam_interaction_imputed_v2 %>% filter(IBD_inflam_interaction_imputed_v2$pval_adj_bh<0.05)
##Extract significant SNP files----
write.table(IBD_inflam_interaction_p_adj_0.05_imputed_v2$variant_id,file = "SNP_IBD_inflam_interaction_p_adj_0.05_imputed.csv",sep="\t",col.names = F, row.names = F,quote = F)

snp_annotation_eqtl<-read.delim("/Users/jingjing/Downloads/Annotation/RUNS/RESULTS_92/snp_annotation.txt")
## 27 eQTL in inflammed tissues.
library(biomaRt)

# Set up the connection to Ensembl-----
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID<-gsub('\\..+$', '', IBD_inflam_interaction_p_adj_0.05_imputed_v2$phenotype_id) ##remove the version id after ".", other wise you can not link those things
IBD_inflam_interaction_p_adj_0.05_imputed_v2$symbol<-mapIds(org.Hs.eg.db,
                   keys = IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID,
                   column = "SYMBOL",
                   keytype = "ENSEMBL")
IBD_inflam_interaction_p_adj_0.05_imputed_v2$GO<-mapIds(org.Hs.eg.db,
                                                            keys = IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID,
                                                            column = "GO",
                                                            keytype = "ENSEMBL")
IBD_inflam_interaction_p_adj_0.05_imputed_v2$ONTOLOGYALL<-mapIds(org.Hs.eg.db,
                                                        keys = IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID,
                                                        column = "ONTOLOGYALL",
                                                        keytype = "ENSEMBL")

IBD_inflam_interaction_p_adj_0.05_imputed_v2$Gene_name<-mapIds(org.Hs.eg.db,
       keys = IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID,
       column = "GENENAME",
       keytype = "ENSEMBL")

library(clusterProfiler)

GO_pathway_Inflammed_eGENE_list<-IBD_inflam_interaction_p_adj_0.05_imputed_v2$b_gi
names(GO_pathway_Inflammed_eGENE_list)<-IBD_inflam_interaction_p_adj_0.05_imputed_v2$ensembl_ID
GO_pathway_Inflammed_eGENE<-gseGO(gene=sort(GO_pathway_Inflammed_eGENE_list,decreasing = TRUE),
      OrgDb=org.Hs.eg.db,
      keyType='ENSEMBL',
      ont="ALL",
      minGSSize = 2, 
      maxGSSize = 800, 
      pvalueCutoff = 0.05, 
      verbose = TRUE,
      pAdjustMethod = "none")
enrich_pathway_Inflammed_eGENE<-enrichGO(gene=names(GO_pathway_Inflammed_eGENE_list),
                                             OrgDb=org.Hs.eg.db,
                                             keyType='ENSEMBL',
                                             ont="ALL",
                                             minGSSize = 2, 
                                             maxGSSize = 800, 
                                             pvalueCutoff = 0.05,
                                             qvalueCutoff = 0.05,
                                             pAdjustMethod = "BH")

View(enrich_pathway_Inflammed_eGENE)
library(SNPannotator)
db <- "1000GENOMES:phase_3:EUR"
server<- "https://rest.ensembl.org"
Inflamed_eSNP_ann<- annotate(IBD_inflam_interaction_p_adj_0.05_imputed_v2$variant_id,
                             server,
                             db,
                             'Inflamed_eSNP_ann_PROXY.xlsx', 
                             LDlist = TRUE, 
                             r2 = 0.5,
                             cadd = FALSE, 
                             geneNames.file = '/Users/jingjing/Documents/R_files/UMCG/Data_base/Gene_Names_Ensembl_104_GRCh38.rds',
                             regulatoryType.file = 
                               '/Users/jingjing/Documents/R_files/UMCG/Data_base/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.rds')


####Extract the eSNP-eGene pairs----
## 1. Host defenses/cell proliferation and differentiation: ENSG00000105697,ENSG00000136869, ENSG00000171049,ENSG00000106799, ENSG00000142871
## 2. Cancer related genes: ENSG00000169035, ENSG00000204614, ENSG00000105926,ENSG00000088882,ENSG00000165071
## 3. Hormone activity/enzyme/receptor for short chain free fatty acids: ENSG00000157017,ENSG00000096395,ENSG00000136881,ENSG00000126262, ENSG00000154646,ENSG0000018475

##1.1Extract pair ENSG00000105697 ----
##extract the expression data for ENSG00000105697 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000105697<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000105697.9") %>% t()
expression_ENSG00000105697 <- as.data.frame(expression_ENSG00000105697)## make a vertical dataframe

colnames(expression_ENSG00000105697)<-expression_ENSG00000105697[4,]
expression_ENSG00000105697$ID <- row.names(expression_ENSG00000105697)

expression_ENSG00000105697<-expression_ENSG00000105697[-c(1:4),]
expression_ENSG00000105697<-expression_ENSG00000105697[,c(2,1)]

##extract the genotype data- rs12610509 from IBD_inflammation samples
rs12610509<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs12610509.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs12610509)
##Merge the genotype data with the expression data
ENSG00000105697_rs12610509<- expression_ENSG00000105697 %>% left_join(rs12610509[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000105697_rs12610509$genotype <- ifelse(ENSG00000105697_rs12610509$A1=="C"& ENSG00000105697_rs12610509$A2=="C",2,
                                                ifelse(ENSG00000105697_rs12610509$A1=="T"& ENSG00000105697_rs12610509$A2=="T",0,1)) ## 0==TT,1=TC,2=CC
View(ENSG00000105697_rs12610509)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000105697_rs12610509 <- ENSG00000105697_rs12610509 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))


##make a scatter plot---- 
ENSG00000105697_rs12610509$ENSG00000105697.9<-as.numeric(ENSG00000105697_rs12610509$ENSG00000105697.9)
library(ggpubr)

P_ENSG00000105697_rs12610509<-
  ggplot(ENSG00000105697_rs12610509, 
         aes(x = genotype, 
             y = log10(ENSG00000105697.9+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | HAMP | rs12610509 | FDR interaction=0.043",
       x = "rs12610509_C",
       y = "Log transformed expression-HAMP",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs12610509_C", 
                     breaks = c(0, 1, 2),
                     label = c("TT", "TC", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000105697_rs12610509.pdf")
P_ENSG00000105697_rs12610509
dev.off()




























##1.2Extract pair ENSG00000136869----
##extract the expression data for ENSG00000136869 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000136869<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000136869.16") %>% t()
expression_ENSG00000136869 <- as.data.frame(expression_ENSG00000136869)## make a vertical dataframe

colnames(expression_ENSG00000136869)<-expression_ENSG00000136869[4,]
expression_ENSG00000136869$ID <- row.names(expression_ENSG00000136869)

expression_ENSG00000136869<-expression_ENSG00000136869[-c(1:4),]
expression_ENSG00000136869<-expression_ENSG00000136869[,c(2,1)]

##Extract the genotype data- rs17397828 from IBD_inflammation samples
rs17397828<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs17397828.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs17397828)
##Merge the genotype data with the expression data
ENSG00000136869_rs17397828<- expression_ENSG00000136869 %>% left_join(rs17397828[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000136869_rs17397828$genotype <- ifelse(ENSG00000136869_rs17397828$A1=="G"& ENSG00000136869_rs17397828$A2=="G",2,
                                              ifelse(ENSG00000136869_rs17397828$A1=="A"& ENSG00000136869_rs17397828$A2=="A",0,1)) ## 0==AA,1=AG,2=GG
View(ENSG00000136869_rs17397828)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000136869_rs17397828 <- ENSG00000136869_rs17397828 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000136869_rs17397828$ENSG00000136869.16<-as.numeric(ENSG00000136869_rs17397828$ENSG00000136869.16)
library(ggpubr)

P_ENSG00000136869_rs17397828<-
  ggplot(ENSG00000136869_rs17397828, 
         aes(x = genotype, 
             y = log(ENSG00000136869.16+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | TLR4 | rs17397828 | FDR interaction=0.018",
       x = "rs17397828_G",
       y = "Log transformed expression-TLR4",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs17397828_G", 
                     breaks = c(0, 1, 2),
                     label = c("AA", "AG", "GG"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000136869_rs17397828.pdf")
P_ENSG00000136869_rs17397828
dev.off()


##1.3Extract pair ENSG00000171049----
##extract the expression data for ENSG00000171049 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000171049<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000171049.9") %>% t()
expression_ENSG00000171049 <- as.data.frame(expression_ENSG00000171049)## make a vertical dataframe

colnames(expression_ENSG00000171049)<-expression_ENSG00000171049[4,]
expression_ENSG00000171049$ID <- row.names(expression_ENSG00000171049)

expression_ENSG00000171049<-expression_ENSG00000171049[-c(1:4),]
expression_ENSG00000171049<-expression_ENSG00000171049[,c(2,1)]

##Extract the genotype data-chr19_51510414_G_GC from IBD_inflammation samples
chr19_51510414_G_GC<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_19:51510414:G:GC.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(chr19_51510414_G_GC)
##Merge the genotype data with the expression data
ENSG00000171049_chr19_51510414_G_GC<- expression_ENSG00000171049 %>% left_join(chr19_51510414_G_GC[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000171049_chr19_51510414_G_GC$genotype <- ifelse(ENSG00000171049_chr19_51510414_G_GC$A1=="GC"& ENSG00000171049_chr19_51510414_G_GC$A2=="GC",2,
                                              ifelse(ENSG00000171049_chr19_51510414_G_GC$A1=="G"& ENSG00000171049_chr19_51510414_G_GC$A2=="G",0,1)) ## 0==GG,1=G/GC,2=GC/GC
View(ENSG00000171049_chr19_51510414_G_GC)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000171049_chr19_51510414_G_GC <- ENSG00000171049_chr19_51510414_G_GC %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000171049_chr19_51510414_G_GC$ENSG00000171049.9<-as.numeric(ENSG00000171049_chr19_51510414_G_GC$ENSG00000171049.9)
library(ggpubr)

P_ENSG00000171049_chr19_51510414_G_GC<-
  ggplot(ENSG00000171049_chr19_51510414_G_GC, 
         aes(x = genotype, 
             y = log(ENSG00000171049.9+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | FPR2 | Chr19_51510414_G_GC | FDR interaction=0.018",
       x = "Chr19_51510414_G_GC",
       y = "Log transformed expression-FPR2",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="Chr19_51510414_G_GC", 
                     breaks = c(0, 1, 2),
                     label = c("GG", "het", "alt"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000171049_chr19_51510414_G_GC.pdf")
P_ENSG00000171049_chr19_51510414_G_GC
dev.off()






##1.4Extract pair ENSG00000106799----
##extract the expression data for ENSG00000106799 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000106799<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000106799.15") %>% t()
expression_ENSG00000106799 <- as.data.frame(expression_ENSG00000106799)## make a vertical data-frame

colnames(expression_ENSG00000106799)<-expression_ENSG00000106799[4,]
expression_ENSG00000106799$ID <- row.names(expression_ENSG00000106799)

expression_ENSG00000106799<-expression_ENSG00000106799[-c(1:4),]
expression_ENSG00000106799<-expression_ENSG00000106799[,c(2,1)]

view(expression_ENSG00000106799)
##Extract the genotype data-rs4407948 from IBD_inflammation samples
rs4407948<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs4407948.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs4407948)
##Merge the genotype data with the expression data
ENSG00000106799_rs4407948<- expression_ENSG00000106799 %>% left_join(rs4407948[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000106799_rs4407948$genotype <- ifelse(ENSG00000106799_rs4407948$A1=="C"& ENSG00000106799_rs4407948$A2=="C",2,
                                                       ifelse(ENSG00000106799_rs4407948$A1=="A"& ENSG00000106799_rs4407948$A2=="A",0,1)) ## 0==AA,1=AC,2=CC
View(ENSG00000106799_rs4407948)

##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000106799_rs4407948 <- ENSG00000106799_rs4407948 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000106799_rs4407948$ENSG00000106799.15<-as.numeric(ENSG00000106799_rs4407948$ENSG00000106799.15)
library(ggpubr)

P_ENSG00000106799_rs4407948<-
  ggplot(ENSG00000106799_rs4407948, 
         aes(x = genotype, 
             y = log(ENSG00000106799.15+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | TGFBR1 | rs4407948 | FDR interaction=0.026",
       x = "rs4407948",
       y = "Log transformed expression-TGFBR1",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs4407948_C", 
                     breaks = c(0, 1, 2),
                     label = c("AA", "AC", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000106799_rs4407948.pdf")
P_ENSG00000106799_rs4407948
dev.off()







##1.5Extract pair ENSG00000142871----
##extract the expression data for ENSG00000142871 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000142871<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000142871.18") %>% t()
expression_ENSG00000142871 <- as.data.frame(expression_ENSG00000142871)## make a vertical data-frame

colnames(expression_ENSG00000142871)<-expression_ENSG00000142871[4,]
expression_ENSG00000142871$ID <- row.names(expression_ENSG00000142871)

expression_ENSG00000142871<-expression_ENSG00000142871[-c(1:4),]
expression_ENSG00000142871<-expression_ENSG00000142871[,c(2,1)]

view(expression_ENSG00000142871)
##Extract the genotype data-rs10782540 from IBD_inflammation samples
rs10782540<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs10782540.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs10782540)
##Merge the genotype data with the expression data
ENSG00000142871_rs10782540<- expression_ENSG00000142871 %>% left_join(rs10782540[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000142871_rs10782540$genotype <- ifelse(ENSG00000142871_rs10782540$A1=="T"& ENSG00000142871_rs10782540$A2=="T",2,
                                             ifelse(ENSG00000142871_rs10782540$A1=="A"& ENSG00000142871_rs10782540$A2=="A",0,1)) ## 0==AA,1=AT,2=TT
View(ENSG00000142871_rs10782540)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000142871_rs10782540 <- ENSG00000142871_rs10782540 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000142871_rs10782540$ENSG00000142871.18<-as.numeric(ENSG00000142871_rs10782540$ENSG00000142871.18)
library(ggpubr)

P_ENSG00000142871_rs10782540<-
  ggplot(ENSG00000142871_rs10782540, 
         aes(x = genotype, 
             y = log(ENSG00000142871.18+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | CCN1 | rs10782540 | FDR interaction=0.009",
       x = "rs10782540",
       y = "Log transformed expression-CCN1",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs10782540_T", 
                     breaks = c(0, 1, 2),
                     label = c("AA", "AT", "TT"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000142871_rs10782540.pdf")
P_ENSG00000142871_rs10782540
dev.off()






##2.1Extract ENSG00000169035(KLK7 have a bit too much missing values)----
##extract the expression data for ENSG00000169035 from the file "filtered_data_qtl_bed_IBD_inflammation"

expression_ENSG00000169035<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000169035.12") %>% t()
expression_ENSG00000169035 <- as.data.frame(expression_ENSG00000169035)## make a vertical data-frame

colnames(expression_ENSG00000169035)<-expression_ENSG00000169035[4,]
expression_ENSG00000169035$ID <- row.names(expression_ENSG00000169035)

expression_ENSG00000169035<-expression_ENSG00000169035[-c(1:4),]
expression_ENSG00000169035<-expression_ENSG00000169035[,c(2,1)]

view(expression_ENSG00000169035)
##Extract the genotype data-chr19_51510414_G_GC from IBD_inflammation samples
chr19_51510414_G_GC<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_19:51510414:G:GC.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(chr19_51510414_G_GC)
##Merge the genotype data with the expression data
ENSG00000169035_chr19_51510414_G_GC<- expression_ENSG00000169035 %>% left_join(chr19_51510414_G_GC[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000169035_chr19_51510414_G_GC$genotype <- ifelse(ENSG00000169035_chr19_51510414_G_GC$A1=="GC"& ENSG00000169035_chr19_51510414_G_GC$A2=="GC",2,
                                              ifelse(ENSG00000169035_chr19_51510414_G_GC$A1=="G"& ENSG00000169035_chr19_51510414_G_GC$A2=="G",0,1)) ## 0==GG,1=G/GC,2=TGC/GC
View(ENSG00000169035_chr19_51510414_G_GC)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000169035_chr19_51510414_G_GC <- ENSG00000169035_chr19_51510414_G_GC %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000169035_chr19_51510414_G_GC$ENSG00000169035.12<-as.numeric(ENSG00000169035_chr19_51510414_G_GC$ENSG00000169035.12)
library(ggpubr)

P_ENSG00000169035_chr19_51510414_G_GC<-
  ggplot(ENSG00000169035_chr19_51510414_G_GC, 
         aes(x = genotype, 
             y = log(ENSG00000169035.12+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | KLK7 | Chr19_51510414_G_GC | FDR interaction=0.002",
       x = "Chr19_51510414_G_GC",
       y = "Log transformed expression-KLK7",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="Chr19_51510414_G_GC", 
                     breaks = c(0, 1, 2),
                     label = c("GG", "het", "alt"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000169035_chr19_51510414_G_GC.pdf")
P_ENSG00000169035_chr19_51510414_G_GC
dev.off()






##2.2Extract ENSG00000204614----
##extract the expression data for ENSG00000204614 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000204614<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000204614.9") %>% t()
expression_ENSG00000204614 <- as.data.frame(expression_ENSG00000204614)## make a vertical data-frame

colnames(expression_ENSG00000204614)<-expression_ENSG00000204614[4,]
expression_ENSG00000204614$ID <- row.names(expression_ENSG00000204614)

expression_ENSG00000204614<-expression_ENSG00000204614[-c(1:4),]
expression_ENSG00000204614<-expression_ENSG00000204614[,c(2,1)]

view(expression_ENSG00000204614)

##Extract the genotype data-rs9278949 from IBD_inflammation samples
rs9278949<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs9278949.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs9278949)
##Merge the genotype data with the expression data
ENSG00000204614_rs9278949<- expression_ENSG00000204614 %>% left_join(rs9278949[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000204614_rs9278949$genotype <- ifelse(ENSG00000204614_rs9278949$A1=="A"& ENSG00000204614_rs9278949$A2=="A",2,
                                                       ifelse(ENSG00000204614_rs9278949$A1=="AAAAAAGAAG"& ENSG00000204614_rs9278949$A2=="AAAAAAGAAG",0,1)) ## 0==ref,1=het,2=AA
View(ENSG00000204614_rs9278949)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000204614_rs9278949 <- ENSG00000204614_rs9278949 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000204614_rs9278949$ENSG00000204614.9<-as.numeric(ENSG00000204614_rs9278949$ENSG00000204614.9)
library(ggpubr)

P_ENSG00000204614_rs9278949<-
  ggplot(ENSG00000204614_rs9278949, 
         aes(x = genotype, 
             y = log(ENSG00000204614.9+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | TRIM40 | rs9278949 | FDR interaction=0.006",
       x = "rs9278949",
       y = "Log transformed expression-TRIM40",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs9278949_A", 
                     breaks = c(0, 1, 2),
                     label = c("ref", "het", "AA"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000204614_rs9278949.pdf")
P_ENSG00000204614_rs9278949
dev.off()









##2.3Extract ENSG00000105926----
##extract the expression data for ENSG00000105926 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000105926<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000105926.16") %>% t()
expression_ENSG00000105926 <- as.data.frame(expression_ENSG00000105926)## make a vertical data-frame

colnames(expression_ENSG00000105926)<-expression_ENSG00000105926[4,]
expression_ENSG00000105926$ID <- row.names(expression_ENSG00000105926)

expression_ENSG00000105926<-expression_ENSG00000105926[-c(1:4),]
expression_ENSG00000105926<-expression_ENSG00000105926[,c(2,1)]

view(expression_ENSG00000105926)

##Extract the genotype data-rs6970988 from IBD_inflammation samples
rs6970988<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs6970988.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs6970988)
##Merge the genotype data with the expression data
ENSG00000105926_rs6970988<- expression_ENSG00000105926 %>% left_join(rs6970988[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000105926_rs6970988$genotype <- ifelse(ENSG00000105926_rs6970988$A1=="A"& ENSG00000105926_rs6970988$A2=="A",2,
                                             ifelse(ENSG00000105926_rs6970988$A1=="C"& ENSG00000105926_rs6970988$A2=="C",0,1)) ## 0==CC,1=CA,2=AA
View(ENSG00000105926_rs6970988)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000105926_rs6970988 <- ENSG00000105926_rs6970988 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000105926_rs6970988$ENSG00000105926.16<-as.numeric(ENSG00000105926_rs6970988$ENSG00000105926.16)
library(ggpubr)

P_ENSG00000105926_rs6970988<-
  ggplot(ENSG00000105926_rs6970988, 
         aes(x = genotype, 
             y = log(ENSG00000105926.16+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | PALS2 | rs6970988 | FDR interaction=0.008",
       x = "rs6970988",
       y = "Log transformed expression-PALS2",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs6970988_A", 
                     breaks = c(0, 1, 2),
                     label = c("CC", "CA", "AA"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000105926_rs6970988.pdf")
P_ENSG00000105926_rs6970988
dev.off()
















##2.4Extract ENSG00000088882----
##extract the expression data for ENSG00000088882 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000088882<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000088882.8") %>% t()
expression_ENSG00000088882 <- as.data.frame(expression_ENSG00000088882)## make a vertical data-frame

colnames(expression_ENSG00000088882)<-expression_ENSG00000088882[4,]
expression_ENSG00000088882$ID <- row.names(expression_ENSG00000088882)

expression_ENSG00000088882<-expression_ENSG00000088882[-c(1:4),]
expression_ENSG00000088882<-expression_ENSG00000088882[,c(2,1)]

view(expression_ENSG00000088882)

##Extract the genotype data-rs79760497 from IBD_inflammation samples
rs79760497<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs79760497.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs79760497)
##Merge the genotype data with the expression data
ENSG00000088882_rs79760497<- expression_ENSG00000088882 %>% left_join(rs79760497[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000088882_rs79760497$genotype <- ifelse(ENSG00000088882_rs79760497$A1=="C"& ENSG00000088882_rs79760497$A2=="C",2,
                                             ifelse(ENSG00000088882_rs79760497$A1=="T"& ENSG00000088882_rs79760497$A2=="T",0,1)) ## 0==TT,1=TC,2=CC
View(ENSG00000088882_rs79760497)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000088882_rs79760497 <- ENSG00000088882_rs79760497 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000088882_rs79760497$ENSG00000088882.8<-as.numeric(ENSG00000088882_rs79760497$ENSG00000088882.8)
library(ggpubr)

P_ENSG00000088882_rs79760497<-
  ggplot(ENSG00000088882_rs79760497, 
         aes(x = genotype, 
             y = log(ENSG00000088882.8+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | CPXM1 | rs79760497 | FDR interaction=0.037",
       x = "rs79760497",
       y = "Log transformed expression-CPXM1",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs79760497_C", 
                     breaks = c(0, 1, 2),
                     label = c("TT", "TC", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000088882_rs79760497.pdf")
P_ENSG00000088882_rs79760497
dev.off()



















##2.5Extract ENSG00000165071----
##extract the expression data for ENSG00000165071 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000165071<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000165071.15") %>% t()
expression_ENSG00000165071 <- as.data.frame(expression_ENSG00000165071)## make a vertical data-frame

colnames(expression_ENSG00000165071)<-expression_ENSG00000165071[4,]
expression_ENSG00000165071$ID <- row.names(expression_ENSG00000165071)

expression_ENSG00000165071<-expression_ENSG00000165071[-c(1:4),]
expression_ENSG00000165071<-expression_ENSG00000165071[,c(2,1)]

view(expression_ENSG00000165071)

##Extract the genotype data-rs2922459from IBD_inflammation samples
rs2922459<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs2922459.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs2922459)
##Merge the genotype data with the expression data
ENSG00000165071_rs2922459<- expression_ENSG00000165071 %>% left_join(rs2922459[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000165071_rs2922459$genotype <- ifelse(ENSG00000165071_rs2922459$A1=="T"& ENSG00000165071_rs2922459$A2=="T",2,
                                              ifelse(ENSG00000165071_rs2922459$A1=="C"& ENSG00000165071_rs2922459$A2=="C",0,1)) ## 0==CC,1=CT,2=TT
View(ENSG00000165071_rs2922459)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000165071_rs2922459 <- ENSG00000165071_rs2922459 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000165071_rs2922459$ENSG00000165071.15<-as.numeric(ENSG00000165071_rs2922459$ENSG00000165071.15)
library(ggpubr)

P_ENSG00000165071_rs2922459<-
  ggplot(ENSG00000165071_rs2922459, 
         aes(x = genotype, 
             y = log(ENSG00000165071.15+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | TMEM71 | rs2922459 | FDR interaction=0.0003",
       x = "rs2922459",
       y = "Log transformed expression-TMEM71",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs2922459_T", 
                     breaks = c(0, 1, 2),
                     label = c("CC", "CT", "TT"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000165071_rs2922459.pdf")
P_ENSG00000165071_rs2922459
dev.off()




























##3.1Extract ENSG00000157017(CHRL missing value a little bit high)----
##extract the expression data for ENSG00000157017 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000157017<- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000157017.16") %>% t()
expression_ENSG00000157017 <- as.data.frame(expression_ENSG00000157017)## make a vertical data-frame

colnames(expression_ENSG00000157017)<-expression_ENSG00000157017[4,]
expression_ENSG00000157017$ID <- row.names(expression_ENSG00000157017)

expression_ENSG00000157017<-expression_ENSG00000157017[-c(1:4),]
expression_ENSG00000157017<-expression_ENSG00000157017[,c(2,1)]

view(expression_ENSG00000157017)

##Extract the genotype data-rs113242288 from IBD_inflammation samples
rs113242288<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs113242288.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs113242288)
##Merge the genotype data with the expression data
ENSG00000157017_rs113242288<- expression_ENSG00000157017 %>% left_join(rs113242288[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000157017_rs113242288$genotype <- ifelse(ENSG00000157017_rs113242288$A1=="C"& ENSG00000157017_rs113242288$A2=="C",2,
                                             ifelse(ENSG00000157017_rs113242288$A1=="T"& ENSG00000157017_rs113242288$A2=="T",0,1)) ## 0==TT,1=TC,2=CC
View(ENSG00000157017_rs113242288)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000157017_rs113242288 <- ENSG00000157017_rs113242288 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000157017_rs113242288$ENSG00000157017.16<-as.numeric(ENSG00000157017_rs113242288$ENSG00000157017.16)
library(ggpubr)

P_ENSG00000157017_rs113242288<-
  ggplot(ENSG00000157017_rs113242288, 
         aes(x = genotype, 
             y = log(ENSG00000157017.16+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | CHRL | rs113242288 | FDR interaction=0.0002",
       x = "rs113242288",
       y = "Log transformed expression-CHRL",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs113242288_C", 
                     breaks = c(0, 1, 2),
                     label = c("TT", "TC", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000157017_rs113242288.pdf")
P_ENSG00000157017_rs113242288
dev.off()





































##3.2Extract ENSG00000096395----
##extract the expression data for ENSG00000096395 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000096395 <- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000096395.11") %>% t()
expression_ENSG00000096395 <- as.data.frame(expression_ENSG00000096395)## make a vertical data-frame

colnames(expression_ENSG00000096395)<-expression_ENSG00000096395[4,]
expression_ENSG00000096395$ID <- row.names(expression_ENSG00000096395)

expression_ENSG00000096395<-expression_ENSG00000096395[-c(1:4),]
expression_ENSG00000096395<-expression_ENSG00000096395[,c(2,1)]

view(expression_ENSG00000096395)

##Extract the genotype data-rs17215168 from IBD_inflammation samples
rs17215168<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs17215168.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs17215168)
##Merge the genotype data with the expression data
ENSG00000096395_rs17215168<- expression_ENSG00000096395 %>% left_join(rs17215168[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000096395_rs17215168$genotype <- ifelse(ENSG00000096395_rs17215168$A1=="A"& ENSG00000096395_rs17215168$A2=="A",2,
                                               ifelse(ENSG00000096395_rs17215168$A1=="T"& ENSG00000096395_rs17215168$A2=="T",0,1)) ## 0==TT,1=TA,2=AA
View(ENSG00000096395_rs17215168)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000096395_rs17215168 <- ENSG00000096395_rs17215168 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000096395_rs17215168$ENSG00000096395.11<-as.numeric(ENSG00000096395_rs17215168$ENSG00000096395.11)
library(ggpubr)

P_ENSG00000096395_rs17215168<-
  ggplot(ENSG00000096395_rs17215168, 
         aes(x = genotype, 
             y = log(ENSG00000096395.11+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | MLN | rs17215168 | FDR interaction=0.002",
       x = "rs17215168",
       y = "Log transformed expression-MLN",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs17215168_A", 
                     breaks = c(0, 1, 2),
                     label = c("TT", "TA", "AA"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000096395_rs17215168.pdf")
P_ENSG00000096395_rs17215168
dev.off()










































##3.3Extract ENSG00000136881.12----
##extract the expression data for ENSG00000136881 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000136881 <- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000136881.12") %>% t()
expression_ENSG00000136881 <- as.data.frame(expression_ENSG00000136881)## make a vertical data-frame

colnames(expression_ENSG00000136881)<-expression_ENSG00000136881[4,]
expression_ENSG00000136881$ID <- row.names(expression_ENSG00000136881)

expression_ENSG00000136881<-expression_ENSG00000136881[-c(1:4),]
expression_ENSG00000136881<-expression_ENSG00000136881[,c(2,1)]

view(expression_ENSG00000136881)

##Extract the genotype data-rs2067538 from IBD_inflammation samples
rs2067538<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs2067538.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs2067538)
##Merge the genotype data with the expression data
ENSG00000136881_rs2067538<- expression_ENSG00000136881 %>% left_join(rs2067538[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000136881_rs2067538$genotype <- ifelse(ENSG00000136881_rs2067538$A1=="A"& ENSG00000136881_rs2067538$A2=="A",2,
                                              ifelse(ENSG00000136881_rs2067538$A1=="C"& ENSG00000136881_rs2067538$A2=="C",0,1)) ## 0==CC,1=CA,2=AA
View(ENSG00000136881_rs2067538)
##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000136881_rs2067538 <- ENSG00000136881_rs2067538 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000136881_rs2067538$ENSG00000136881.12<-as.numeric(ENSG00000136881_rs2067538$ENSG00000136881.12)
library(ggpubr)

P_ENSG00000136881_rs2067538<-
  ggplot(ENSG00000136881_rs2067538, 
         aes(x = genotype, 
             y = log(ENSG00000136881.12+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | BAAT | rs2067538 | FDR interaction=0.0003",
       x = "rs2067538",
       y = "Log transformed expression-BAAT",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs2067538_A", 
                     breaks = c(0, 1, 2),
                     label = c("CC", "CA", "AA"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000136881_rs2067538.pdf")
P_ENSG00000136881_rs2067538
dev.off()















































##3.4Extract ENSG00000126262.5----
##extract the expression data for ENSG00000126262 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000126262 <- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000126262.5") %>% t()
expression_ENSG00000126262 <- as.data.frame(expression_ENSG00000126262)## make a vertical data-frame

colnames(expression_ENSG00000126262)<-expression_ENSG00000126262[4,]
expression_ENSG00000126262$ID <- row.names(expression_ENSG00000126262)

expression_ENSG00000126262<-expression_ENSG00000126262[-c(1:4),]
expression_ENSG00000126262<-expression_ENSG00000126262[,c(2,1)]

view(expression_ENSG00000126262)

##Extract the genotype data-rs12610509 from IBD_inflammation samples
rs12610509<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_rs12610509.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(rs12610509)

##Merge the genotype data with the expression data
ENSG00000126262_rs12610509<- expression_ENSG00000126262 %>% left_join(rs12610509[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000126262_rs12610509$genotype <- ifelse(ENSG00000126262_rs12610509$A1=="C"& ENSG00000126262_rs12610509$A2=="C",2,
                                             ifelse(ENSG00000126262_rs12610509$A1=="T"& ENSG00000126262_rs12610509$A2=="T",0,1)) ## 0==TT,1=TC,2=CC
View(ENSG00000126262_rs12610509)

##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000126262_rs12610509 <- ENSG00000126262_rs12610509 %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000126262_rs12610509$ENSG00000126262.5<-as.numeric(ENSG00000126262_rs12610509$ENSG00000126262.5)
library(ggpubr)

P_ENSG00000126262_rs12610509<-
  ggplot(ENSG00000126262_rs12610509, 
         aes(x = genotype, 
             y = log(ENSG00000126262.5+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | FFAR2 | rs12610509 | FDR interaction=0.026",
       x = "rs12610509",
       y = "Log transformed expression-FFAR2",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="rs12610509_C", 
                     breaks = c(0, 1, 2),
                     label = c("TT", "TC", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000126262_rs12610509.pdf")
P_ENSG00000126262_rs12610509
dev.off()




















































##3.5Extract ENSG00000154646.5----
##extract the expression data for ENSG00000154646 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000154646 <- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000154646.9") %>% t()
expression_ENSG00000154646 <- as.data.frame(expression_ENSG00000154646)## make a vertical data-frame

colnames(expression_ENSG00000154646)<-expression_ENSG00000154646[4,]
expression_ENSG00000154646$ID <- row.names(expression_ENSG00000154646)

expression_ENSG00000154646<-expression_ENSG00000154646[-c(1:4),]
expression_ENSG00000154646<-expression_ENSG00000154646[,c(2,1)]

view(expression_ENSG00000154646)

##Extract the genotype data-Chr21_18264601_CAAAA_C from IBD_inflammation samples
Chr21_18264601_CAAAA_C<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_21:18264601:CAAAA:C.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(Chr21_18264601_CAAAA_C)

##Merge the genotype data with the expression data
ENSG00000154646_Chr21_18264601_CAAAA_C<- expression_ENSG00000154646 %>% left_join(Chr21_18264601_CAAAA_C[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000154646_Chr21_18264601_CAAAA_C$genotype <- ifelse(ENSG00000154646_Chr21_18264601_CAAAA_C$A1=="C"& ENSG00000154646_Chr21_18264601_CAAAA_C$A2=="C",2,
                                              ifelse(ENSG00000154646_Chr21_18264601_CAAAA_C$A1=="CAAAA"& ENSG00000154646_Chr21_18264601_CAAAA_C$A2=="CAAAA",0,1)) ## 0==ref,1=het,2=CC
View(ENSG00000154646_Chr21_18264601_CAAAA_C)

##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000154646_Chr21_18264601_CAAAA_C <- ENSG00000154646_Chr21_18264601_CAAAA_C %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000154646_Chr21_18264601_CAAAA_C$ENSG00000154646.9<-as.numeric(ENSG00000154646_Chr21_18264601_CAAAA_C$ENSG00000154646.9)
library(ggpubr)
View(ENSG00000154646_Chr21_18264601_CAAAA_C)

P_ENSG00000154646_Chr21_18264601_CAAAA_C<-
  ggplot(ENSG00000154646_Chr21_18264601_CAAAA_C, 
         aes(x = genotype, 
             y = log(ENSG00000154646.9+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | TMPRSS15 | Chr21_18264601_CAAAA_C | FDR interaction=0.002",
       x = "Chr21_18264601_CAAAA_C",
       y = "Log transformed expression-TMPRSS15",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="Chr21_18264601_CAAAA_C", 
                     breaks = c(0, 1, 2),
                     label = c("ref", "het", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000154646_Chr21_18264601_CAAAA_C.pdf")
P_ENSG00000154646_Chr21_18264601_CAAAA_C
dev.off()






















































##3.6Extract ENSG00000184752----
##extract the expression data for ENSG00000184752 from the file "filtered_data_qtl_bed_IBD_inflammation"
expression_ENSG00000184752 <- filtered_data_qtl_bed_IBD_inflammation %>% filter(phenotype_id == "ENSG00000184752.14") %>% t()
expression_ENSG00000184752 <- as.data.frame(expression_ENSG00000184752)## make a vertical data-frame

colnames(expression_ENSG00000184752)<-expression_ENSG00000184752[4,]
expression_ENSG00000184752$ID <- row.names(expression_ENSG00000184752)

expression_ENSG00000184752<-expression_ENSG00000184752[-c(1:4),]
expression_ENSG00000184752<-expression_ENSG00000184752[,c(2,1)]

view(expression_ENSG00000184752)

##Extract the genotype data-Chr12_94630256_CA_C from IBD_inflammation samples
Chr12_94630256_CA_C<-read_delim("/Users/jingjing/Documents/R_files/RNA-seq/eQTL_imputed/eSNP_eGene_pair_v2/IBD_inflam_imputed_12:94630256:CA:C.ped",col_names = c("IID","FID","PID","MID","SEX","AGE","A1","A2"))
View(Chr12_94630256_CA_C)

##Merge the genotype data with the expression data
ENSG00000184752_Chr12_94630256_CA_C<- expression_ENSG00000184752 %>% left_join(Chr12_94630256_CA_C[,c(2,7,8)],by=c("ID"="FID"))
ENSG00000184752_Chr12_94630256_CA_C$genotype <- ifelse(ENSG00000184752_Chr12_94630256_CA_C$A1=="C"& ENSG00000184752_Chr12_94630256_CA_C$A2=="C",2,
                                                          ifelse(ENSG00000184752_Chr12_94630256_CA_C$A1=="CA"& ENSG00000184752_Chr12_94630256_CA_C$A2=="CA",0,1)) ## 0==ref,1=het,2=CC
View(ENSG00000184752_Chr12_94630256_CA_C)

##Merge with the smoking data----
View(Inflam_IBD_qtl) ##meta data
ENSG00000184752_Chr12_94630256_CA_C <- ENSG00000184752_Chr12_94630256_CA_C %>% left_join(Inflam_IBD_qtl[,c(1,22)],by=c("ID"="#IID"))
##Make a scatter plot---- 
ENSG00000184752_Chr12_94630256_CA_C$ENSG00000184752.14<-as.numeric(ENSG00000184752_Chr12_94630256_CA_C$ENSG00000184752.14)
library(ggpubr)
View(ENSG00000184752_Chr12_94630256_CA_C)

P_ENSG00000184752_Chr12_94630256_CA_C<-
  ggplot(ENSG00000184752_Chr12_94630256_CA_C, 
         aes(x = genotype, 
             y = log(ENSG00000184752.14+1),
             group=Smoking_status,
             color = Smoking_status,)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = TRUE,aes(fill=Smoking_status),alpha=0.2,show.legend = F)+
  labs(title = "IBD inflammation | NDUFA12 | Chr12_94630256_CA_C | FDR interaction=0.006",
       x = "Chr12_94630256_CA_C",
       y = "Log transformed expression-NDUFA12",
       color = "Smoking status") +
  theme_minimal() +
  scale_color_manual(values = c("Current" = "#DC0000FF", "Never" = "#3C5488FF","Former"="#00A087FF"))+
  scale_x_continuous(name="Chr12_94630256_CA_C", 
                     breaks = c(0, 1, 2),
                     label = c("ref", "het", "CC"))+
  theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold",margin = margin(0,0,15,0)),
        plot.margin = unit(c(1,0,0,0),"cm"),
        legend.position = "bottom")

pdf("P_ENSG00000184752_Chr12_94630256_CA_C.pdf")
P_ENSG00000184752_Chr12_94630256_CA_C
dev.off()






























































##combine----
ggarrange_ieqtl_inflam_IBD<-ggarrange(P_ENSG00000105697_rs12610509,
                                      P_ENSG00000136869_rs17397828,
                                      P_ENSG00000171049_chr19_51510414_G_GC,
                                      P_ENSG00000106799_rs4407948,
                                      P_ENSG00000142871_rs10782540,
                                      P_ENSG00000169035_chr19_51510414_G_GC,
                                      P_ENSG00000204614_rs9278949,
                                      P_ENSG00000105926_rs6970988,
                                      P_ENSG00000088882_rs79760497,
                                      P_ENSG00000165071_rs2922459,
                                      P_ENSG00000157017_rs113242288,
                                      P_ENSG00000096395_rs17215168,
                                      P_ENSG00000136881_rs2067538,
                                      P_ENSG00000126262_rs12610509,
                                      P_ENSG00000154646_Chr21_18264601_CAAAA_C,
                                      P_ENSG00000184752_Chr12_94630256_CA_C,
                                      ncol=4,nrow=4,
                                      labels = c("A","B","C","D",
                                                 "E","F","G","H",
                                                 "I","J","K","L",
                                                 "M","N","O","P"))

ggsave("i_eqtl_inflam_IBD_IMPUTED.pdf", ggarrange_ieqtl_inflam_IBD, width = 30, height = 30, dpi = 600,limitsize = FALSE)

dev.off()




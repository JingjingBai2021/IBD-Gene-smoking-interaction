
##########################################################
##########Prepare the input files  for tensorqtl##########
##########################################################
# Load GTF and extract gene coordinates
Smoking_RNA_seq_closest_ADD_merged_meta<-readxl::read_xlsx("Smoking_RNA_seq_closest_ADD_merged_meta.xlsx")
bed<-as.data.frame(txi_qc[["abundance"]]) ##set TPM as the the phenotype_id
bed$gene_id<-rownames(bed)
bed<-bed[,c(731,1:730)]

gtf <- import("/Users/jingjing/ClusterDrives/gearshift_groups/umcg-weersma/tmp01/Jingjing/tools/Salmon/gencode.v46.annotation.gtf")
gene_coords <- as.data.frame(gtf[gtf$type == "gene"]) %>%
  dplyr::select( chr = seqnames, gene_id = gene_id, start, end) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(start = end - 1)

merged_bed_IBD <- bed%>%
  left_join(gene_coords, by = "gene_id")

##Merged_bed_IBD file contains the RNA-seq for all the IBD patients who have smoking data.
write.table(
  merged_bed_IBD,
  file = "merged_bed_IBD.bed",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = FALSE)

view(merged_bed_IBD)

######tensorqtl only allow one genotype to one expression data.frame, so have to select only one sample from each patient to make the bed file###
##Split the bed file based on inflammation
###Inflammation-Bed-ibd file

Inflam_IBD<-Smoking_RNA_seq_closest_ADD_merged_meta %>% filter(Inflammation=="Yes") 
Inflam_IBD_qtl<-  as.data.frame(unique(Inflam_IBD$`Research ID`))
colnames(Inflam_IBD_qtl)<-"Research ID" 
Inflam_IBD_qtl<-Inflam_IBD_qtl%>% left_join(Inflam_IBD,by="Research ID",multiple = "first")
####Prepare the input files for non-inflammation IBD eqtl analysis

Non_inflam_IBD<-Smoking_RNA_seq_closest_ADD_merged_meta %>% filter(Inflammation=="No") 
Non_inflam_IBD_qtl<- as.data.frame(unique(Non_inflam_IBD$`Research ID`))
colnames(Non_inflam_IBD_qtl)<-"Research ID" 
Non_inflam_IBD_qtl<-Non_inflam_IBD_qtl%>% left_join(Non_inflam_IBD,by="Research ID",multiple = "first")

##link the genotype ID with those RNA-seq ID 
UMCG_Smoking<-readxl::read_xlsx("UMCG_Smoking.xlsx",sheet = "Sheet4")

Inflam_IBD_qtl<- Inflam_IBD_qtl %>% left_join(UMCG_Smoking[,c(1,2)], by="Research ID")
Inflam_IBD_qtl<-Inflam_IBD_qtl[,c(58,1:57)]
Inflam_IBD_qtl <- Inflam_IBD_qtl %>% filter(!is.na(Inflam_IBD_qtl$`#IID`)) ###180 samples with genotype data and smoking and rna-seq:)

Non_inflam_IBD_qtl<-Non_inflam_IBD_qtl %>% left_join(UMCG_Smoking[,c(1,2)], by="Research ID")
Non_inflam_IBD_qtl<-Non_inflam_IBD_qtl[,c(58,1:57)]
Non_inflam_IBD_qtl<-Non_inflam_IBD_qtl %>% filter(!is.na(Non_inflam_IBD_qtl$`#IID`)) ## 265 samples with genotype data and smoking and rna-seq data

write.table(Non_inflam_IBD_qtl[,c(1,1)],file="Non_inflam_IBD_id.txt",quote=F,row.names = F,col.names = T,sep="\t") ##write the table to subset the genotype data 
##plink2 --bfile EUR --make-bed --keep Non_inflam_IBD_id.txt --output-chr chrM --out IBD_Non_Inflam_38_chr

##subset the bed file--Inflammation
qtl_bed_IBD_inflammation <- merged_bed_IBD %>% select("chr","start","end","gene_id",Inflam_IBD_qtl$RNA_seq_ID)
colnames(qtl_bed_IBD_inflammation)[5:184] <- Inflam_IBD_qtl$`#IID`
colnames(qtl_bed_IBD_inflammation)[4]<-"phenotype_id"

##subset the bed file--Non-Inflammation
qtl_bed_IBD_Non_inflammation <- merged_bed_IBD %>% select("chr","start","end","gene_id",Non_inflam_IBD_qtl$RNA_seq_ID)
colnames(qtl_bed_IBD_Non_inflammation)[5:269] <- Non_inflam_IBD_qtl$`#IID` ## convert the RNA-seq ID to genotype IDs  

##Remove the genes with low expression levels (TPM < 0.1 in more than 80% samples)

low_tpm_count_inflam <- rowSums(filtered_data_qtl_bed_IBD_inflammation[,c(5:184)] < 0.1)
filtered_data_qtl_bed_IBD_inflammation <- qtl_bed_IBD_inflammation %>%
  filter(low_tpm_count <= 0.8 * 180)

low_tpm_count_non_inflam <- rowSums(qtl_bed_IBD_Non_inflammation[,c(5:269)] < 0.1)
filtered_data_qtl_bed_IBD_Non_inflammation <- qtl_bed_IBD_Non_inflammation %>%
  filter(low_tpm_count_non_inflam <= 0.8 * 265)


## read genotype fam file 
IBD_inflam_fam<-read_delim("IBD_inflam.text",col_names = c("FID","IID","A","B","C","D")) ##168 samples
IBD_non_inflam_fam<-read_delim("IBD_Non_Inflam_38_chr.txt",col_names = c("FID","IID","A","B","C","D"))## 259 samples

filtered_data_qtl_bed_IBD_Non_inflammation <- filtered_data_qtl_bed_IBD_Non_inflammation %>% select("chr","start","end","gene_id",IBD_non_inflam_fam$FID)
filtered_data_qtl_bed_IBD_inflammation <- filtered_data_qtl_bed_IBD_inflammation %>% select("chr","start","end","phenotype_id",IBD_inflam_fam$FID)

colnames(filtered_data_qtl_bed_IBD_Non_inflammation)[4]<-"phenotype_id"

write.table(
  filtered_data_qtl_bed_IBD_inflammation,
  file = "filtered_data_qtl_bed_IBD_inflammation.bed",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = FALSE
) ## write expression bed file for inflammed samples

write.table(
  filtered_data_qtl_bed_IBD_Non_inflammation,
  file = "filtered_data_qtl_bed_IBD_Non_inflammation.bed",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = FALSE
) ## write expression bed file for non-inflamed samples


colnames(filtered_data_qtl_bed_IBD_inflammation)[5:184]
write.table(
  colnames(filtered_data_qtl_bed_IBD_inflammation)[5:184],
  file = "qtl_bed_IBD_inflammation.csv",
  sep = "\t",
  quote = FALSE,
  col.names = F,
  row.names = FALSE)                                         

##covariate file need to be dummy coded for categorical variables
covariate<-Inflam_IBD_qtl[,c(1,19,22,23,26,27,28,36)] ## extract the possible confounders

covariate$Smoking_status <- ifelse(covariate$Smoking_status =="Never",0,
                                   ifelse(covariate$Smoking_status =="Current",1,2))
covariate$Batch <- ifelse(covariate$Batch =="2016",1,2)
covariate$Location_rough <- ifelse(covariate$Location_rough =="colon",1,2)
covariate$sex<- as.numeric(covariate$sex)
covariate$Diagnosis<- ifelse(covariate$Diagnosis =="CD",1,2)
covariate$biological_use<-as.numeric(covariate$biological_use)
covariate<-covariate[order(covariate$`#IID`),]
which(is.na(covariate$biological_use))
covariate[137,8]<-0 ##missing

covariate<-t(covariate)
colnames(covariate)<-covariate[1,]
covariate<-covariate[-1,]
covariate<- as.data.frame(covariate)
covariate <-covariate %>% select(IBD_inflam_fam$FID)

write.table(
  covariate,
  file = "covariate_Inflam_IBD2.csv",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = T
)                                              

write.table(
  t(covariate[2,]),
  file = "Interaction_Inflam_IBD.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = T
)   

View(covariate_Non_Inflam)
covariate_Non_Inflam<-Non_inflam_IBD_qtl[,c(1,19,22,23,26,27,28,36),] ## extract the possible confounders, add the biological use
covariate_Non_Inflam$Smoking_status <- ifelse(covariate_Non_Inflam$Smoking_status =="Never",0,
                                   ifelse(covariate_Non_Inflam$Smoking_status =="Current",1,2))
covariate_Non_Inflam$Batch <- ifelse(covariate_Non_Inflam$Batch =="2016",1,2)
covariate_Non_Inflam$Location_rough <- ifelse(covariate_Non_Inflam$Location_rough =="colon",1,2)
covariate_Non_Inflam$sex<- as.numeric(covariate_Non_Inflam$sex)
covariate_Non_Inflam$Diagnosis<- ifelse(covariate_Non_Inflam$Diagnosis =="CD",1,2)
covariate_Non_Inflam$biological_use<-as.numeric(covariate_Non_Inflam$biological_use)
covariate_Non_Inflam<-covariate_Non_Inflam[order(covariate_Non_Inflam$`#IID`),]
which(is.na(covariate_Non_Inflam[,8]))
covariate_Non_Inflam[189,8]<-0


covariate_Non_Inflam<-t(covariate_Non_Inflam)
colnames(covariate_Non_Inflam)<-covariate_Non_Inflam[1,]
covariate_Non_Inflam<-covariate_Non_Inflam[-1,]
covariate_Non_Inflam<- as.data.frame(covariate_Non_Inflam)
covariate_Non_Inflam <-covariate_Non_Inflam %>% select(IBD_non_inflam_fam$FID)

non_inflam_ibd_id<-colnames(covariate_Non_Inflam) %>% as.data.frame()
colnames(non_inflam_ibd_id)[1]<-"#FID"
write.table(
  non_inflam_ibd_id,
  file = "non_inflam_ibd_id.csv",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = F
)   




write.table(
  covariate_Non_Inflam,
  file = "covariate_Non_Inflam_IBD2.csv",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = T
)   

write.table(
  t(covariate_Non_Inflam[2,]),
  file = "Interaction_Non_Inflam_IBD.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = T,
  row.names = T
) ##prepare the interaction term file








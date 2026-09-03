
##Current_never-----
setwd("/Users/jingjing/Documents/R_files/Hetero_CD_UC_swap/Current_CD_UC/")
##Hetero_CD_UC_swap folder is the next version of opposite variants by swapping the major allele as effect allele
##Extract the variants both are significant in CD an UC (P<0.01)----
library(readr)
library(metafor)
Meta_current_CD<- readr::read_table("CD.meta")
Meta_current_UC<-readr::read_table("UC.meta")
Overlapping_current<- Reduce(intersect,
                             list(Meta_current_CD[Meta_current_CD$P<0.01,]$SNP,
                                  Meta_current_UC[Meta_current_UC$P<0.01,]$SNP))
Overlapping_current <- as.data.frame(Overlapping_current)
colnames(Overlapping_current)[1]<- "SNP"

##Add CD effect size----
library(tidyverse)
Overlapping_current <- Overlapping_current %>% 
  left_join(Meta_current_CD[,c(1:4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_current)[5]<- "CD_P"
colnames(Overlapping_current)[6]<- "CD_OR"
colnames(Overlapping_current)[7]<- "CD_Q"

##Add UC effect size----
Overlapping_current <- Overlapping_current %>% left_join(Meta_current_UC[,c(3,4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_current)[9]<- "UC_P"
colnames(Overlapping_current)[10]<- "UC_OR"

##Extract the overlapped variants that have opposite effect in CD and UC----
Opposite_current<-Overlapping_current[Overlapping_current$CD_OR>1.0 & Overlapping_current$UC_OR< 1.0,]
Opposite_current2<-Overlapping_current[Overlapping_current$CD_OR<1.0 & Overlapping_current$UC_OR> 1.0,]
Opposite_current<-rbind(Opposite_current,Opposite_current2)

##Add position and ALT allele and 95% CI from summary statistics----
##CD
Meta_current_CD$SE <- abs((log(Meta_current_CD$OR))/qnorm(Meta_current_CD$P/2))
Meta_current_CD$"Lower" <- round(exp(log(Meta_current_CD$OR)-1.96*Meta_current_CD$SE),digits = 2)
Meta_current_CD$"Upper" <- round(exp(log(Meta_current_CD$OR)+1.96*Meta_current_CD$SE),digits = 2)
Meta_current_CD$"95% CI" <- paste(format(Meta_current_CD$"Lower",digits=2),format(Meta_current_CD$"Upper",digits=2), sep = "-")
Meta_current_CD$"OR (95% CI)" <- paste0(format(Meta_current_CD$OR,digits = 2)," (", Meta_current_CD$`95% CI`, ")")
Meta_current_CD$"OR (95% CI)"<- str_replace_all(Meta_current_CD$"OR (95% CI)"," ","")
Meta_current_CD$"95% CI" <-NULL

##UC
Meta_current_UC$SE <- abs((log(Meta_current_UC$OR))/qnorm(Meta_current_UC$P/2))
Meta_current_UC$"Lower" <- round(exp(log(Meta_current_UC$OR)-1.96*Meta_current_UC$SE),digits = 2)
Meta_current_UC$"Upper" <- round(exp(log(Meta_current_UC$OR)+1.96*Meta_current_UC$SE),digits = 2)
Meta_current_UC$"95% CI" <- paste(format(Meta_current_UC$"Lower",digits=2),format(Meta_current_UC$"Upper",digits=2), sep = "-")
Meta_current_UC$"OR (95% CI)" <- paste0(format(Meta_current_UC$OR,digits = 2)," (", Meta_current_UC$`95% CI`, ")")
Meta_current_UC$"OR (95% CI)"<- str_replace_all(Meta_current_UC$"OR (95% CI)"," ","")
Meta_current_UC$"95% CI" <-NULL

##add 95%CI to opposite variants
Opposite_current<- Opposite_current %>% left_join(Meta_current_CD[,c(3,15,16,17,18)], by= c("SNP"="SNP"))
colnames(Opposite_current)
Opposite_current<-Opposite_current[c(1,2,3,4,5,6,7,12,13,14,15,9,10,11)]

colnames(Opposite_current)[8]<-"CD_SE"
colnames(Opposite_current)[9] <- "CD_Lower"
colnames(Opposite_current)[10] <- "CD_Upper"
colnames(Opposite_current)[11] <- "CD_OR(95%CI)"

Opposite_current<- Opposite_current %>% left_join(Meta_current_UC[,c(3,15,16,17,18)], by= c("SNP"="SNP"))

colnames(Opposite_current)[14] <- "UC_Q"
colnames(Opposite_current)[15] <- "UC_SE"
colnames(Opposite_current)[16] <- "UC_Lower"
colnames(Opposite_current)[17] <- "UC_Upper"
colnames(Opposite_current)[18] <- "UC_OR(95%CI)"
colnames(Opposite_current)[4]<-"ALT"

Opposite_current$`UC_OR(95%CI)`<- paste0(format(Opposite_current$UC_OR,digits = 2)," (", paste(format(Opposite_current$UC_Lower,digits=2),format(Opposite_current$UC_Upper,digits=2), sep = "-"), ")")
Opposite_current$`UC_OR(95%CI)`<- str_replace_all(Opposite_current$`UC_OR(95%CI)`," ","")


##remove the missing value lines
Opposite_current <- drop_na(Opposite_current)

##Forest plot all contrast----

Opposite_current$ratio_CD_UC<-(Opposite_current$CD_OR)/(Opposite_current$UC_OR)
Opposite_current_filter<- Opposite_current %>% filter(ratio_CD_UC>10|ratio_CD_UC<0.1)

##Annotation
library(topr)
Opposite_current_filter<-annotate_with_nearest_gene(
  Opposite_current_filter %>% dplyr::rename(CHROM=CHR,POS=BP),
  protein_coding_only = FALSE,
  build = 38,
  .chr_map = NULL
)


##Former_never----
setwd("/Users/jingjing/Documents/R_files/Hetero_CD_UC_swap/Former_CD_UC/")
##Hetero_CD_UC_swap folder is the next version of opposite variants by swapping the major allele as effect allele
##Extract the variants both are significant in CD an UC (P<0.01)----
library(readr)
Meta_former_CD<- readr::read_table("CD.meta")
Meta_former_UC<-readr::read_table("UC.meta")
Overlapping_former<- Reduce(intersect,
                            list(Meta_former_CD[Meta_former_CD$P<0.01,]$SNP,
                                 Meta_former_UC[Meta_former_UC$P<0.01,]$SNP))
Overlapping_former <- as.data.frame(Overlapping_former)
colnames(Overlapping_former)[1]<- "SNP"

##Add CD effect size----
library(tidyverse)
Overlapping_former <- Overlapping_former %>% 
  left_join(Meta_former_CD[,c(1:4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_former)[5]<- "CD_P"
colnames(Overlapping_former)[6]<- "CD_OR"
colnames(Overlapping_former)[7]<- "CD_Q"

##Add UC effect size----
Overlapping_former <- Overlapping_former %>% left_join(Meta_former_UC[,c(3,4,7,9,11)], by = c("SNP"="SNP"))
colnames(Overlapping_former)[9]<- "UC_P"
colnames(Overlapping_former)[10]<- "UC_OR"

##Extract the overlapped variants that have opposite effect in CD and UC----
Opposite_former<-Overlapping_former[Overlapping_former$CD_OR>1.0 & Overlapping_former$UC_OR< 1.0,]
Opposite_former2<-Overlapping_former[Overlapping_former$CD_OR<1.0 & Overlapping_former$UC_OR> 1.0,]
Opposite_former<-rbind(Opposite_former,Opposite_former2)

##Add position and ALT allele and 95% CI from summary statistics----
##CD
Meta_former_CD$SE <- abs((log(Meta_former_CD$OR))/qnorm(Meta_former_CD$P/2))
Meta_former_CD$"Lower" <- round(exp(log(Meta_former_CD$OR)-1.96*Meta_former_CD$SE),digits = 2)
Meta_former_CD$"Upper" <- round(exp(log(Meta_former_CD$OR)+1.96*Meta_former_CD$SE),digits = 2)
Meta_former_CD$"95% CI" <- paste(format(Meta_former_CD$"Lower",digits=2),format(Meta_former_CD$"Upper",digits=2), sep = "-")
Meta_former_CD$"OR (95% CI)" <- paste0(format(Meta_former_CD$OR,digits = 2)," (", Meta_former_CD$`95% CI`, ")")
Meta_former_CD$"OR (95% CI)"<- str_replace_all(Meta_former_CD$"OR (95% CI)"," ","")
Meta_former_CD$"95% CI" <-NULL

##UC
Meta_former_UC$SE <- abs((log(Meta_former_UC$OR))/qnorm(Meta_former_UC$P/2))
Meta_former_UC$"Lower" <- round(exp(log(Meta_former_UC$OR)-1.96*Meta_former_UC$SE),digits = 2)
Meta_former_UC$"Upper" <- round(exp(log(Meta_former_UC$OR)+1.96*Meta_former_UC$SE),digits = 2)
Meta_former_UC$"95% CI" <- paste(format(Meta_former_UC$"Lower",digits=2),format(Meta_former_UC$"Upper",digits=2), sep = "-")
Meta_former_UC$"OR (95% CI)" <- paste0(format(Meta_former_UC$OR,digits = 2)," (", Meta_former_UC$`95% CI`, ")")
Meta_former_UC$"OR (95% CI)"<- str_replace_all(Meta_former_UC$"OR (95% CI)"," ","")
Meta_former_UC$"95% CI" <-NULL

##add 95%CI to opposite variants
Opposite_former<- Opposite_former %>% left_join(Meta_former_CD[,c(3,16,17,18)], by= c("SNP"="SNP"))
colnames(Opposite_former)
Opposite_former<-Opposite_former[c(1,2,3,4,5,6,7,12,13,14,9,10,11)]
colnames(Opposite_former)[8] <- "CD_Lower"
colnames(Opposite_former)[9] <- "CD_Upper"
colnames(Opposite_former)[10] <- "CD_OR(95%CI)"

Opposite_former<- Opposite_former %>% left_join(Meta_former_UC[,c(3,16,17,18)], by= c("SNP"="SNP"))
colnames(Opposite_former)[14] <- "UC_Lower"
colnames(Opposite_former)[15] <- "UC_Upper"
colnames(Opposite_former)[16] <- "UC_OR(95%CI)"
colnames(Opposite_former)[4]<-"ALT"

##remove the missing value lines
Opposite_former <- drop_na(Opposite_former)

Opposite_former$`UC_OR(95%CI)`<- paste0(format(Opposite_former$UC_OR,digits = 2)," (", paste(format(Opposite_former$UC_Lower,digits=2),format(Opposite_former$UC_Upper,digits=2), sep = "-"), ")")
Opposite_former$`UC_OR(95%CI)`<- str_replace_all(Opposite_former$`UC_OR(95%CI)`," ","")

Opposite_former$`CD_OR(95%CI)`<- paste0(format(Opposite_former$CD_OR,digits = 2)," (", paste(format(Opposite_former$CD_Lower,digits=2),format(Opposite_former$CD_Upper,digits=2), sep = "-"), ")")

Opposite_former$`CD_OR(95%CI)`<- str_replace_all(Opposite_former$`CD_OR(95%CI)`," ","")

Opposite_former$ratio_CD_UC<-(Opposite_former$CD_OR)/(Opposite_former$UC_OR)

View(Opposite_former_filter)
Opposite_former_filter<- Opposite_former %>% filter(ratio_CD_UC>10|ratio_CD_UC<0.1)

Opposite_former_filter<-annotate_with_nearest_gene(
  Opposite_former_filter %>% dplyr::rename(CHROM=CHR,POS=BP),
  protein_coding_only = FALSE,
  build = 38,
  .chr_map = NULL
)



##Ever_never----
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

Opposite_ever$`UC_OR(95%CI)`<- paste0(format(Opposite_ever$UC_OR,digits = 2)," (", paste(format(Opposite_ever$UC_Lower,digits=2),format(Opposite_ever$UC_Upper,digits=2), sep = "-"), ")")
Opposite_ever$`UC_OR(95%CI)`<- str_replace_all(Opposite_ever$`UC_OR(95%CI)`," ","")

##remove the missing value lines
Opposite_ever <- drop_na(Opposite_ever)

Opposite_ever$ratio_CD_UC<-(Opposite_ever$CD_OR)/(Opposite_ever$UC_OR)
Opposite_ever_filter<- Opposite_ever %>% filter(ratio_CD_UC>10|ratio_CD_UC<0.1)
View(Opposite_ever_filter)

Opposite_ever_filter<-annotate_with_nearest_gene(
  Opposite_ever_filter %>% dplyr::rename(CHROM=CHR,POS=BP),
  protein_coding_only = FALSE,
  build = 38,
  .chr_map = NULL)



##combine all contrasts
Combine_opposite<- rbind(arrange(Opposite_ever_filter,CHROM,POS),#n=8
                         arrange(Opposite_current_filter,CHROM,POS),#n=12
                         arrange(Opposite_former_filter,CHROM,POS))#n=10


colnames(Combine_opposite)[17]<-"Effect size ratio"

Combine_opposite$"Chr:Pos" <- paste0(Combine_opposite$CHROM,":",Combine_opposite$POS)

Combine_opposite<-Combine_opposite[c(20,1,18,4,
                                     10,5,6,8,9,
                                     16,11,12,14,15,17)]
View(Combine_opposite)


Combine_opposite[nrow(Combine_opposite)+1,]<- NA ##ADD new rows, three times
Combine_opposite<-Combine_opposite[c(31,1:8,
                                     32,9:20,
                                     33,21:30),]

Combine_opposite[1,]$`Chr:Pos`<- paste0("Ever/Never")
Combine_opposite[10,]$`Chr:Pos`<- paste0("Current/Never")
Combine_opposite[23,]$`Chr:Pos`<- paste0("Former/Never")

Combine_opposite$CD_P<- format(Combine_opposite$CD_P,scientific = T)
Combine_opposite$UC_P<-format(Combine_opposite$UC_P,scientific = T)


Combine_opposite$`Effect size ratio`<-round(Combine_opposite$`Effect size ratio`,digits = 2)

Combine_opposite$`Effect size ratio` <- as.character(Combine_opposite$`Effect size ratio`)
Combine_opposite[c(1,10,23),c(2:6,10,11,15)] <- ""


Combine_opposite$""<- paste(rep(" ", 30), collapse = " ")

colnames(Combine_opposite)[4]<-"Major allele"
colnames(Combine_opposite)[5]<-"OR(95%CI).CD"
colnames(Combine_opposite)[6]<-"P.CD"
colnames(Combine_opposite)[10]<-"OR(95%CI).UC"
colnames(Combine_opposite)[11]<-"P.UC"
colnames(Combine_opposite)[3]<-"Nearest gene"


##
library(forestploter)
library(grid)

tm<- forest_theme(base_size = 10,
                  refline_gp = gpar(lyt="solid"),
                  ci_pch = c(15,18),
                  ci_col = c("#023e8a", "#D21312"),
                  legend_name = "Group",
                  legend_value = c("CD", "UC"),
                  legend_position = "right",
                  vertline_lty = c("dotted"),
                  vertline_col = c("#bababa"))

png('Opposite_extremes_forest_swap_new.png', res = 600, width = 15, height = 15, units = "in")


forest(Combine_opposite[,c(1, 2, 3, 4, 16, 5,6,10,11)],
       est = list(Combine_opposite$CD_OR,
                  Combine_opposite$UC_OR),
       lower = list(Combine_opposite$CD_Lower,
                    Combine_opposite$UC_Lower), 
       upper = list(Combine_opposite$CD_Upper,
                    Combine_opposite$UC_Upper),
       ci_column = c(5),
       ref_line = 1,
       clip = c(0.001, 5),
       vert_line = c(0.5,1.5),
       xlim = c(0,5),
       nudge_y = 0.2,
       theme = tm)
dev.off()

write.table(Combine_opposite,"Combine_opposite.csv")
##Add heterogeneity test estimates----


lines <- readLines("Combine_opposite.csv")
lines <- gsub("\r", "", lines)           # remove Windows carriage return
lines <- gsub('""', '"', lines)          # collapse double-quotes to single
lines <- gsub('^"|"$', '', lines)        # strip leading/trailing quote per line

# Write cleaned lines to temp file and read as space-separated
tmp <- tempfile()
writeLines(lines, tmp)

df <- read.table(tmp, 
                 header = TRUE, 
                 stringsAsFactors = FALSE, 
                 check.names = FALSE,
                 quote = '"',
                 sep = " ",
                 fill = TRUE)

# Drop the empty trailing column
df <- df[, 1:15]

colnames(df)
View(df)

write.csv(df, "Combine_opposite2.csv", row.names = FALSE, quote = FALSE)




##Add estimates for the heterogeneity tests.
results <- data.frame(
  SNP = df$SNP,
  Q   = NA,
  p_Q = NA,
  I2  = NA
)

for (i in 1:nrow(df)) {
  
  log_OR1 <- log(df$CD_OR[i])
  log_OR2 <- log(df$UC_OR[i])
  
  SE1 <- (log(df$CD_Upper[i]) - log(df$CD_Lower[i])) / (2 * 1.96)
  SE2 <- (log(df$UC_Upper[i]) - log(df$UC_Lower[i])) / (2 * 1.96)
  
  w1 <- 1 / SE1^2
  w2 <- 1 / SE2^2
  
  mean_log_OR <- (w1 * log_OR1 + w2 * log_OR2) / (w1 + w2)
  
  Q  <- w1 * (log_OR1 - mean_log_OR)^2 + w2 * (log_OR2 - mean_log_OR)^2
  
  results$Q[i]   <- round(Q, 2)
  results$p_Q[i] <- formatC(pchisq(Q, df = 1, lower.tail = FALSE),format = "e",digits=3)
  results$I2[i]  <- round(max(0, (Q - 1) / Q) * 100, 2)
}

print(results)



##Add statistics
df_new_1<-df[1:9,] %>% left_join(results[1:9,],by="SNP")
print(df_new_1)
df_new_2<-df[10:22,] %>% left_join(results[10:22,],by="SNP")
print(df_new_2)
df_new_3<-df[23:33,] %>% left_join(results[23:33,],by="SNP")
print(df_new_3)

df_new<-rbind(df_new_1,df_new_2,df_new_3)

write.csv(df_new,file = "Opposite_SNPs_CD_UC_with_Q.csv",row.names = FALSE, quote = FALSE)

colnames(df_new)
df_new[16]<-NULL

##Forestplots----
library(forestploter)
library(grid)
print(df_new)

df_new$P.CD<-formatC(df_new$P.CD, format = "e", digits = 3)
tm<- forest_theme(base_size = 10,
                  refline_gp = gpar(lyt="solid"),
                  ci_pch = c(15,18),
                  ci_col = c("#023e8a", "#D21312"),
                  legend_name = "Group",
                  legend_value = c("CD", "UC"),
                  legend_position = "right",
                  vertline_lty = c("dotted"),
                  vertline_col = c("#bababa"))

df_new$""<- paste(rep(" ", 30), collapse = " ")

colnames(df_new)

p<-forest(df_new[,c(1, 2, 3, 4, 19, 5,6,10,11,16:18)],
       est = list(df_new$CD_OR,
                  df_new$UC_OR),
       lower = list(df_new$CD_Lower,
                    df_new$UC_Lower), 
       upper = list(df_new$CD_Upper,
                    df_new$UC_Upper),
       ci_column = c(5),
       ref_line = 1,
       clip = c(0.001, 5),
       vert_line = c(0.5,1.5),
       xlim = c(0,5),
       nudge_y = 0.2,
       theme = tm)
ggsave("Opposite_snps_Q.pdf", p, width = 18, height = 10, dpi = 600,limitsize = FALSE)




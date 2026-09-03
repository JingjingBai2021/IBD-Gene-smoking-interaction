library(circlize)
library(tidyverse)
##Prepare the input files------
#########Ever/never#############################------
##IBD----
IBD_ever_lead<- Meta_IBD_ever_lead[,c(1,4,3)]
IBD_ever_lead$value <- 1

##CD----
setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_ever")
Meta_CD_ever_lead<-readr::read_table("second_clumped_lead.csv.clumped")
CD_ever_lead<-Meta_CD_ever_lead[,c(1,4,3)]
colnames(CD_ever_lead)[1]<-"chr"
colnames(CD_ever_lead)[2]<-"start"
CD_ever_lead$"end"<-CD_ever_lead$start
CD_ever_lead$chr<-paste("chr",CD_ever_lead$chr,sep = "")
CD_ever_lead$value <-1
CD_ever_lead<- CD_ever_lead[,c(1,2,4,3,5)]

##To make the plot look complete, add some invalid values
CD_ever_lead[96,1]<-"chr20"
CD_ever_lead[97,1]<-"chr22"
CD_ever_lead[98,1]<-"chr22"
CD_ever_lead[96,2]<-0
CD_ever_lead[97,2]<-0
CD_ever_lead[98,2]<-0

##UC----
setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_ever")
Meta_UC_ever_lead<-readr::read_table("second_clumped_lead.csv.clumped")
UC_ever_lead<- Meta_UC_ever_lead[,c(1,4,3)]
colnames(UC_ever_lead)[1]<-"chr"
colnames(UC_ever_lead)[2]<-"start"
UC_ever_lead$"end"<-UC_ever_lead$start
UC_ever_lead$chr<- paste("chr",UC_ever_lead$chr,sep="")
UC_ever_lead$value<-1
UC_ever_lead<-UC_ever_lead[,c(1,2,4,3,5)]

UC_ever_lead[50,1]<-"chr16"
UC_ever_lead[52,1]<-"chr18"
UC_ever_lead[52,1]<-"chr21"
UC_ever_lead[53,1]<-"chr2"

UC_ever_lead[50,2]<-0
UC_ever_lead[51,2]<-0
UC_ever_lead[52,2]<-0
UC_ever_lead[53,2]<-0

UC_ever_lead[50,3]<-0
UC_ever_lead[51,3]<-0
UC_ever_lead[52,3]<-0
UC_ever_lead[53,3]<-0


#############Current/never----
##IBD----
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_current/")
Meta_IBD_current_lead<-readr::read_table("second_clumped_lead.csv.clumped")
IBD_current_lead<-Meta_IBD_current_lead[,c(1,4,3)]
colnames(IBD_current_lead)[1]<-"chr"
colnames(IBD_current_lead)[2]<-"start"
IBD_current_lead$"end"<-IBD_current_lead$start
IBD_current_lead$chr<-paste("chr",IBD_current_lead$chr,sep = "")
IBD_current_lead$value <-1
IBD_current_lead<- IBD_current_lead[,c(1,2,4,3,5)]

nrow(IBD_current_lead)
IBD_current_lead[72,1]<-"chr20"
IBD_current_lead[72,2]<-0
IBD_current_lead[72,3]<-0

##CD----
setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_current/")
Meta_CD_current_lead<-readr::read_table("second_clumped_lead.csv.clumped")
CD_current_lead<-Meta_CD_current_lead[,c(1,4,3)]
colnames(CD_current_lead)[1]<-"chr"
colnames(CD_current_lead)[2]<-"start"
CD_current_lead$"end"<-CD_current_lead$start
CD_current_lead$chr<-paste("chr",CD_current_lead$chr,sep = "")
CD_current_lead$value <-1
CD_current_lead<- CD_current_lead[,c(1,2,4,3,5)]

table(CD_current_lead$chr)
nrow(CD_current_lead)
CD_current_lead[82,1]<-"chr19"
CD_current_lead[83,1]<-"chr20"
CD_current_lead[84,1]<-"chr21"

CD_current_lead[82,2]<-0
CD_current_lead[83,2]<-0
CD_current_lead[84,2]<-0

CD_current_lead[82,3]<-0
CD_current_lead[83,3]<-0
CD_current_lead[84,3]<-0

##UC----
setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_current/")
Meta_UC_current_lead<-readr::read_table("second_clumped_lead.csv.clumped")
UC_current_lead<-Meta_UC_current_lead[,c(1,4,3)]
colnames(UC_current_lead)[1]<-"chr"
colnames(UC_current_lead)[2]<-"start"
UC_current_lead$"end"<-UC_current_lead$start
UC_current_lead$chr<-paste("chr",UC_current_lead$chr,sep = "")
UC_current_lead$value <-1
UC_current_lead<- UC_current_lead[,c(1,2,4,3,5)]

table(UC_current_lead$chr)
nrow(UC_current_lead)

UC_current_lead[87,1]<-"chr18"
UC_current_lead[87,2]<-0
UC_current_lead[87,3]<-0

################Former_never----
##IBD----
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_former/")
Meta_IBD_former_lead<-readr::read_table("second_clumped_lead.csv.clumped")
IBD_former_lead<-Meta_IBD_former_lead[,c(1,4,3)]
colnames(IBD_former_lead)[1]<-"chr"
colnames(IBD_former_lead)[2]<-"start"
IBD_former_lead$"end"<-IBD_former_lead$start
IBD_former_lead$chr<-paste("chr",IBD_former_lead$chr,sep = "")
IBD_former_lead$value <-1
IBD_former_lead<- IBD_former_lead[,c(1,2,4,3,5)]

table(IBD_former_lead$chr)
nrow(IBD_former_lead)
IBD_former_lead[74,1]<-"chr20"
IBD_former_lead[75,1]<-"chr21"

IBD_former_lead[74,2]<-0
IBD_former_lead[75,2]<-0

IBD_former_lead[74,3]<-0
IBD_former_lead[75,3]<-0

##CD----
setwd("/Users/jingjing/Documents/R_files/Meta/META_CD_former/")
Meta_CD_former_lead<-readr::read_table("second_clumped_lead.csv.clumped")
CD_former_lead<-Meta_CD_former_lead[,c(1,4,3)]
colnames(CD_former_lead)[1]<-"chr"
colnames(CD_former_lead)[2]<-"start"
CD_former_lead$"end"<-CD_former_lead$start
CD_former_lead$chr<-paste("chr",CD_former_lead$chr,sep = "")
CD_former_lead$value <-1
CD_former_lead<- CD_former_lead[,c(1,2,4,3,5)]

table(CD_former_lead$chr)
nrow(CD_former_lead)

CD_former_lead[78,1]<-"chr14"
CD_former_lead[79,1]<-"chr15"

CD_former_lead[78,2]<-0
CD_former_lead[79,2]<-0
CD_former_lead[78,3]<-0
CD_former_lead[79,3]<-0
##UC----
setwd("/Users/jingjing/Documents/R_files/Meta/META_UC_former/")
Meta_UC_former_lead<-readr::read_table("second_clumped_lead.csv.clumped")
UC_former_lead<-Meta_UC_former_lead[,c(1,4,3)]
colnames(UC_former_lead)[1]<-"chr"
colnames(UC_former_lead)[2]<-"start"
UC_former_lead$"end"<-UC_former_lead$start
UC_former_lead$chr<-paste("chr",UC_former_lead$chr,sep = "")
UC_former_lead$value <-1
UC_former_lead<- UC_former_lead[,c(1,2,4,3,5)]

table(UC_former_lead$chr)
nrow(UC_former_lead)

UC_former_lead[50,1]<-"chr12"
UC_former_lead[51,1]<-"chr14"
UC_former_lead[52,1]<-"chr19"
UC_former_lead[53,1]<-"chr5"

UC_former_lead[50,2]<-0
UC_former_lead[51,2]<-0
UC_former_lead[52,2]<-0
UC_former_lead[53,2]<-0
UC_former_lead[50,3]<-0
UC_former_lead[51,3]<-0
UC_former_lead[52,3]<-0
UC_former_lead[53,3]<-0

##labels----
setwd("/Users/jingjing/Documents/R_files/Meta/META_IBD_ever/Yadav")
labels<- read_delim("Cir_label.csv")
labels<-labels[,c(2,3,4,1,5,6)]
colnames(labels)[1]<-"chr"
colnames(labels)[2]<-"start"
labels$"end"<-labels$Pos
labels<-labels[,c(1,2,7,3:6)]
labels$chr<-paste("chr",labels$chr,sep="")
library(writexl)
write_xlsx(labels,"labels.xlsx")
##Initialize the plot----
##set parameters
circos.par("start.degree" = 85,
           "gap.degree"=1,
           "cell.padding"= c(0,0,0,0))


##load the hg38 as outside circular
circos.initializeWithIdeogram(species = "hg38",
                              plotType = NULL,
                              chromosome.index= paste0("chr", c(1:22,"Y")))


##add labels outside the track
circos.genomicLabels(labels, 
                     labels.column = 4, 
                     side = "outside",
                     cex=0.15,
                     line_lwd = c(0.1),
                     labels_height=c(0.01),
                     track.margin = c(0,0),
                     padding = c(0.06),
                     col = c('lightblue', '#DBB854')[as.numeric(labels$Group)],
                     line_col = c('lightblue', '#DBB854')[as.numeric(labels$Group)])

##change the chromosome index number format and position
circos.track(ylim = c(0, 1), 
             panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = "white")
  circos.text(mean(xlim), mean(ylim), 
              gsub(".*chr", "", CELL_META$sector.index),
              cex = 0.3, 
              col = "black",
              facing = "inside", 
              niceFacing = TRUE)
}, track.height = 0.04, bg.border = NA,track.margin=c(0.01,0))


##generate a new tract for IBD_ever
circos.genomicTrackPlotRegion(IBD_ever_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0.005),
                              panel.fun = function(region, value, ...) {
  circos.rect(col = "#FBEFE2",border = NA)
  circos.genomicLines(region, value, type = "h",col = "#ca6702")
},bg.border = NA)

###CD_ever
circos.genomicTrackPlotRegion(CD_ever_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#FBEFE2",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#ca6702")
                              },bg.border = NA)


###UC_ever
circos.genomicTrackPlotRegion(UC_ever_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#FBEFE2",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#ca6702")
                              },bg.border = NA)
##IBD_current
circos.genomicTrackPlotRegion(IBD_current_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0.01),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#F8F1F7",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#C57AA2")
                              },bg.border = NA)

##CD_current
circos.genomicTrackPlotRegion(CD_current_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#F8F1F7",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#C57AA2")
                              },bg.border = NA)

##UC_current
circos.genomicTrackPlotRegion(UC_current_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#F8F1F7",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#C57AA2")
                              },bg.border = NA)


##IBD_former
circos.genomicTrackPlotRegion(IBD_former_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0.01),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#EAF5F4",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#6b9080")
                              },bg.border = NA)
##CD_former
circos.genomicTrackPlotRegion(CD_former_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#EAF5F4",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#6b9080")
                              },bg.border = NA)

##UC_former
circos.genomicTrackPlotRegion(UC_former_lead, 
                              ylim = c(0, 1),
                              track.height = 0.04,
                              track.margin =c(0.0015,0),
                              panel.fun = function(region, value, ...) {
                                circos.rect(col = "#EAF5F4",border = NA)
                                circos.genomicLines(region, value, type = "h",col = "#6b9080")
                              },bg.border = NA)

circos.clear()



# Copyright statement comment
# Author comment
    # Chord visualization of 9 sample PBMC FACs data with metabolic and immune markers
# File description comment, including purpose of program, inputs, and outputs
    # Get the chord visualization  comparing the association/mean metabolic marker MFI across all cell types

# Input: FACs summarized output for each sample; Output: Chord graph with legend
# Check if all the required packages are installed
package_list <- c("readxl", "dplyr", "circlize", "ComplexHeatmap")
package_new <- package_list[!(package_list %in% installed.packages()[, "Package"])]
if(length(package_new))
  install.packages(package_new)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!("ComplexHeatmap" %in% installed.packages()[, "Package"]))
  BiocManager::install("ComplexHeatmap")

# library() statements to load the packages required
library(readxl)
library(dplyr)
# chord visualization
library(circlize)
library(ComplexHeatmap)  
####################################################################################################
# set the working directory to source file location where the data is located in the same folder
# setwd("Path to file")
####################################################################################################
# Read the gMFI file
DT <- read_xlsx("gMFI_data_Untreated.xlsx", sheet = "Untreated")
marker_group <- read_xlsx("gMFI_data_Untreated.xlsx", sheet = "Sheet1")
immune <- pull(marker_group[, 1])
immune <- as.vector(na.omit(immune)) 
metabolic <- pull(marker_group[, 2]) # get the names of metabolic markers
metabolic <- as.vector(na.omit(metabolic)) 


# clean the names of the gMFI dataset
markers <- sapply(strsplit(colnames(DT)[3:ncol(DT)], " gMFI "), function(x) x[2])
markers <- gsub("[\\(\\)]", "", markers)
colnames(DT) <- c("Celltype", "Donor", markers)
DT <- as.data.frame(DT[, c(colnames(DT)[1:2], immune, metabolic)])
# Replace negative values with 0
DT[, 3:ncol(DT)] <- apply(DT[, 3:ncol(DT)], 2, function(x) {x[x<0] <- 0;x})
# View(DT)
# Simplify the cell names
DT$Celltype <- as.factor(DT$Celltype)
levels(DT$Celltype) <- c("B.IgD-IgM-", "B.IgD+", "B.IgD+IgM+", "B.IgM+", "mDCs.CD11c+", 
                         "Monocytes.CD16-", "Monocytes.CD16+", "NK.CD56+CD16-", "NK.CD56+CD16+", 
                         "NK.CD56+CD3+", "NK.CD56hiCD16-", "pDCs.CD123+", "T.CD4+", "T.CD8+"  
)

###############################################################################
# Get the spearman correlation between a cell type vs all others ( e.g. B vs all Others) 
# for differential metabolic marker MFI expression

# compute the spearsman correlation 
cor.r <- NULL
for (i in unique(DT$Celltype)){
  temp <- DT[, metabolic]
  temp$cell <- DT$Celltype == i
  cor.r <- rbind(cor.r, c(i, cor(temp, method = "spearman")["cell", ]))
}
row.names <- cor.r[, 1] # get the cell types names
cor.r <- cor.r[, -c(1, ncol(cor.r))] # extract the correlation values only
cor.r<-apply(cor.r, 2, as.numeric)
rownames(cor.r)<-row.names
write.csv(cor.r, "Unactivated Spearman 201906.csv")  ## save the spearman correlation

## Plot the Chord Diagram
col_fun_corr <- colorRamp2(c(-1,  0, 1),  c("#0000FF",  "#FFFFFF", "#FF0000"),  transparency = 0)
lgd_links <- Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_corr, title = "cor")

grid.col <- c("B.IgD-IgM-" = "darkgreen",  "B.IgD+" = "darkgreen",  "B.IgD+IgM+" = "darkgreen", 
             "B.IgM+" = "darkgreen",  "mDCs.CD11c+" = "darkgreen",  "Monocytes.CD16-" = "darkgreen", 
             "Monocytes.CD16+" =  "darkgreen",  "NK.CD56+CD16-" =  "darkgreen",  "NK.CD56+CD16+" = "darkgreen", 
             "NK.CD56+CD3+" = "darkgreen",  "NK.CD56hiCD16-" = "darkgreen",  "pDCs.CD123+" = "darkgreen", 
             "T.CD4+" = "darkgreen", "T.CD8+" = "darkgreen",  
             "GLUT1" = "darkgrey",  "SLC20A1" = "darkgrey",  "ASS1" = "darkgrey",  "IDH2" = "darkgrey" , 
             "G6PD" = "darkgrey",  "ACAC" = "darkgrey",  "PRDX2" = "darkgrey",  "HK1" = "darkgrey",  
             "CPT1A" = "darkgrey",  "ATP5A" = "darkgrey"
)
png(paste("Correlation spearman", "Unactivated.png"), height = 1000, width = 1200, res = 150)

chordDiagram(cor.r,  grid.col = grid.col,  annotationTrack = "grid",  col = col_fun_corr,transparency = 0)
circos.track(track.index = 1,  panel.fun = function(x,  y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
 circos.text(mean(xlim),  ylim[1],  sector.name,  facing = "clockwise",  
             niceFacing = TRUE,  adj = c(0.2,  0.5), cex = 0.8) # labelling
},  bg.border = NA)


draw(lgd_links, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))

dev.off()

##############################################################################
# Read the gMFI file of T cells 
DT <- read_xlsx("AllCD4_Act_Arc 201907.xlsx", sheet = "Sheet1")
marker_group <- read_xlsx("AllCD4_Act_Arc 201907.xlsx", sheet = "Sheet2")
immune <- pull(marker_group[, 1])
immune <- as.vector(na.omit(immune)) 
metabolic <- pull(marker_group[, 2]) # get the names of metabolic markers
metabolic <- as.vector(na.omit(metabolic)) 

# spearman correlation
cor.r<-cor(DT[,-1], method = "spearman")
cor.Immuno.Meta<-cor.r[immune,metabolic]
write.csv(cor.Immuno.Meta,"Spearman correlation of MFI between activated and unactivated.csv")

## Plot the Chord Diagram
grid.col = c("CD69" = "darkgreen",  "CCR7" = "darkgreen",  "HLA-DR" = "darkgreen", 
             "CD25" = "darkgreen",  "CD45RA" = "darkgreen",   
             "ASS1" = "darkgrey",  "IDH2" = "darkgrey" , 
             "G6PD" = "darkgrey",  "ACAC" = "darkgrey",  "PRDX2" = "darkgrey",  "HK1" = "darkgrey",  
             "CPT1A" = "darkgrey", "GLUT1" = "darkgrey",  "ATP5A" = "darkgrey"
)
# continuous Legend
lgd_links = Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_corr, 
                   title = "cor")


png(paste("spearman corr of MFI difference between", "activated vs unactivated legend.png"), height = 1000, width = 1200, res = 150)
chordDiagram(cor.Immuno.Meta,  grid.col = grid.col,  annotationTrack = "grid",  col = col_fun_corr, transparency = 0)
circos.track(track.index = 1,  panel.fun = function(x,  y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim),  ylim[1],  sector.name,  facing = "clockwise",  
               niceFacing = TRUE,  adj = c(0.2,  0.5), cex = 0.5) # labelling
},  bg.border = NA)

draw(lgd_links, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))


dev.off()

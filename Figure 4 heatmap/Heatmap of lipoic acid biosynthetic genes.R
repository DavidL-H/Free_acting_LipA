# This is a script for reading the comparative genome analysis performed in Microbesonline, analyzing it and 
# creating plots


# Load standard libraries
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyr)
library(reshape2)
library(stringi)
library(readxl)
library(dplyr)
library(plyr)
library(Biobase)

dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
# dir_name <- "C:/Your_filepath_here"
setwd(dir_name)

#Read the excel file
# All lipoic acid gene occurences
LAcomgen <- read_excel("LipA_B_lplA_LipM_LipL_gcvH orthologs.xlsx")

# Metadata for only bacterial genomes
bac_meta <- read.csv("MO_bacteria_metadata.csv", header = T)

# Lipoic acid gene occurences and metadata for ONLY bacterial genomes
Alldat <- merge(LAcomgen,bac_meta,by="Genome")

# Subset to only include ones with >0 LipAs in genomes.
Alldat <- subset(Alldat, Alldat$ecLipA > 0)
# Remove un-needed columns for plotting heatmap
heatdata <- Alldat[,-c(3:10,20:38)]

# Make heatmap
datamatrix <- as.matrix(heatdata[,c(3:8,10)])
rownames(datamatrix)= Alldat$Genome

# Create own palatte for coloring the heatmap
my_palette <- colorRampPalette(c("#ff5c5c", "#fff75c", "#5cff64"))(n = 14)


col_breaks = c(seq(0,0.4,length=5),               # for red
               seq(0.5,1.5,length=5),           # for yellow
               seq(2,7,length=5))             # for green



#Plote the heatmap, full A4 paige


#pdf("AllgenomesLAmetabolismHeatmap.pdf",onefile=T, paper='A4', height = 11, width = 8)
svg("AllgenomesLAmetabolismHeatmap.svg",onefile=T, height = 11, width = 8)

heatmap.2(datamatrix, Colv=FALSE, Rowv=TRUE, dendrogram='row',
          main = "Putative Lipoic acid metabolism genes \nin 1062 LipA encoding bacterial genomes", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),   # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,
          scale="none",
          cexRow=0.055,         # Font size row labels
          cexCol=1,             # Font size col labels
          lhei = c(1,8),        #Adjusts the size of the color key
          lwid = c(1,4),        #Adjusts the size of the color key
          srtCol = 0.5,         #Makes the column labels horizontal instead of vertical
          adjCol = c(NA, 0));    #Makes the column label align with columns

# enable color transition at specified limits)

dev.off()
openPDF(normalizePath("AllgenomesLAmetabolismHeatmap.pdf"))

heatmap.2(datamatrix)

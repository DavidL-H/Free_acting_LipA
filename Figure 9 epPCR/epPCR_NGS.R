# 210401
# David Lennox-Hvenekilde
# Analysis of deep-sequenced error-Prone PCR libraries, selected and un-selected for functionality in E. coli
# Genome Sequencer Illumina HiSeq, NovaSeq 6000 S2 PE150 XP. The ecoLipA ORF was sequenced at a coverage of 
# ~4*105. CLC genomics workbench 11.0 was used to analyze data. Reads were mapped with standard settings to 
# the ecoLipA ORF plus ~100 bp up and downstream. Basic variant detection tool was used, setting "Minimum coverage" 
# and "Minimum count" to 1 and "Minimum frequency" to 10-8, with all "Noise filters" disabled.
#################################################################################
# Load libraries and data
library(stringr)
library(plyr)
library(ggplot2)
library(seqinr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(ggsci)
library(plotly)
library(ggpubr)
library(processx)

# Clear the workspace. set the working dir. Make sure that this script is in the same folder as the .csv data files
rm(list = ls(all.names = TRUE))
dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)

# Load the selected and unselected plasmid library variant detection results, the output from CLC genomics.
Sel<-read.csv("Selected epPCR ecoLipA deep sequencing.csv")
NSel<-read.csv("Unselected epPCR ecoLipA deep sequencing.csv")
NSel<-NSel[,-c(16,17,18,19,20,21,22)]
# We performed the plasmid deep sequencing with including ~100 bp on each side of the lipA ORF to improve mapping
# Now we are only intersted in the ORF, which we filter here:
Sel <- subset(Sel, Sel$Overlapping.annotations=="CDS: ecoLipA")
NSel <- subset(NSel, NSel$Overlapping.annotations=="CDS: ecoLipA")

# Due to the very low frequency cut-off for variants detection, there are likely many sequencing artifacts
# This practically means that many of the detected variants in the library are not true variants, and there
# Needs to be some kind of cut-off for the frequency of detection to limit these.
# This will be done by using the frame-shift mutations in the selected library as a kind of internal standard.
# This can be done under the assumption that a frame-shift in the ORF of LipA will always destroy its function,
# And thus a true frame-shift will never be found in the selected LipA epPCR library.
# 
hist(Sel$Count)
hist(NSel$Count)

# Look at the distribution of frame-shift frequencies
# A frame shift is defined as an insertion or deletion of a number of nucleotides not in triplicate codon format
unique(Sel$Type)
Sel_FS <- subset(Sel, (Sel$Type=="Deletion" | Sel$Type=="Insertion"))
Sel_nonFS <- subset(Sel, !(Sel$Type=="Deletion" | Sel$Type=="Insertion"))

# Plot the cumulative ditribution functions to visualize the differences between FS and non FS mutations
plot(ecdf(Sel_nonFS$Count),cex=0, xlim = c(0,1000), col="#4DBBD5", main="", xlab="Count", ylab="Fn(X)")
lines(ecdf(Sel_FS$Count),
      col="#3C5488", cex=0)
abline(h=0.95, col="#E64B35", lty=2, lwd= 0.5)
abline(v=69, col="#E64B35", lty=2, lwd= 0.5)

# Look at the function of the cumulative distribution of the frame
Fn <- ecdf(Sel_FS$Count)
# Find the cutoff for the cumulative function of the frame-shift for the selected dataset
Fn(139)  # returns the percentiles for x
Fn(69)
# Using 139 as the count (coverage) cut-off removes 99% of these false positives frame-shifts and 32% of non-frameshift
# Using 69 as the count (coverage) cut-off removes 95% of these false positives frame-shifts and 10% of non-frameshift
Fn_2 <- ecdf(Sel_nonFS$Count)
Fn_2(139)
Fn_2(69)

# Calucating ratio of mutability##############################################################
# Comparing specific residues
# First we are cleaning up the dataframes a bit more by extracting the mutations into a new column
NSel$AA<-sub('.*\\.', '', NSel$Amino.acid.change)
Sel$AA<-sub('.*\\.', '', Sel$Amino.acid.change)

# Now subset these with the chosen "count" cut-off value and NOT frame-shift
NSel_red<-subset(NSel, (NSel$Count > 69) & !(NSel$Type=="Deletion" | NSel$Type=="Insertion"))
Sel_red<-subset(Sel, (Sel$Count > 69) & !(Sel$Type=="Deletion" | Sel$Type=="Insertion"))

# Remove silent mutations, i.e. no aminoacid change in empty rows
NSel_red<-subset(NSel_red, !(NSel_red$AA==""))
Sel_red<-subset(Sel_red, !(Sel_red$AA==""))

# Extract the amino acid position of each change and put it in a new column
NSel_red$AApos<-extract_numeric(NSel_red$AA)
Sel_red$AApos<-extract_numeric(Sel_red$AA)

# Extract the aminoacid of the ORF that has been mutated into a new column
NSel_red$AAORF<-str_extract(NSel_red$AA, "[aA-zZ]+")
Sel_red$AAORF<-str_extract(Sel_red$AA, "[aA-zZ]+")

# Now we have a quick look at the mutations of each residue in ecLipA
# Which ones are more often mutated in the selected residue vs. the unselected experiment
LipAlen<-322
unique(NSel_red$AApos)


# Make a new DF and populate it with data on residue mutability
Res_mut = data.frame(matrix("", ncol = 6, nrow = LipAlen))
colnames(Res_mut)<- c("position", "res", "Freq_Sel", "Freq_NSel", "Count_Sel", "Count_NSel")
Res_mut$position<-as.integer(row.names(Res_mut))

for (n in 1:nrow(Res_mut)){
  Res_mut$Freq_Sel[n]<-as.integer(nrow(subset(Sel_red, Sel_red$AApos==n)))
  Res_mut$Freq_NSel[n]<-as.integer(nrow(subset(NSel_red, NSel_red$AApos==n)))
  Res_mut$res[n]<-unique(c(subset(Sel_red, Sel_red$AApos==n)$AAORF, subset(NSel_red, NSel_red$AApos==n)$AAORF))
  Res_mut$Count_Sel[n] <- as.integer(sum(subset(Sel_red, Sel_red$AApos==n)$Count))
  Res_mut$Count_NSel[n] <- as.integer(sum(subset(NSel_red, NSel_red$AApos==n)$Count))
}

# For some reason the values are strings, here they are reconverted to integers
Res_mut$Count_Sel<- as.integer(Res_mut$Count_Sel)
Res_mut$Count_NSel<-as.integer(Res_mut$Count_NSel)
Res_mut$Freq_Sel<- as.integer(Res_mut$Freq_Sel)
Res_mut$Freq_NSel<-as.integer(Res_mut$Freq_NSel)
Res_mut$position<-as.integer(Res_mut$position)

# Make a new column with the log2 ration of the difference
Res_mut$Ratio<-log2(Res_mut$Freq_NSel/Res_mut$Freq_Sel)
Res_mut$RatioCount<-log2(Res_mut$Count_NSel/Res_mut$Count_Sel)

# Lets have a look at the ones that have a ratio of Non-selected to selected higher than the
# highest of the other way around. The rational is that there is no reason why a residue should be less
# mutable in the non-selected subset (negative log2 ratio), anything that is negative is the result of
# noise in the data and retracting this noise may give us a better overview of what we are actually 
# looking for, i.e., residues that are signifcantly less mutable in the selected subset.

Res_mut_SIG<-subset(Res_mut,Res_mut$RatioCount > -min(Res_mut$RatioCount))
#write.csv(Res_mut_SIG, "Res_mut_SIG_69CO.csv")
write.csv(Res_mut, "epPCR_Res_mut.csv")

# Plot the data ###################################################################
# Here we plot the log2 ratio of the difference in mutability for all residues in ecLipA:
p1 <- plot_ly(x = Res_mut$position, y = Res_mut$RatioCount,  type = 'scatter', fill = 'tozeroy',  mode = "lines", 
              fillcolor = '#8491B4',
              line = list(
                color = '#3C5488'
              ))%>%
  
  add_annotations(x = Res_mut_SIG$position,
                  y = Res_mut_SIG$RatioCount,
                  text = paste(Res_mut_SIG$res, Res_mut_SIG$position, sep=""),
                  xref = "x",
                  yref = "y",
                  showarrow = TRUE,
                  arrowhead = 4,
                  arrowsize = .5,
                  ax = 20,
                  ay = -20,
                  font = list(color = '#E64B35',
                              family = 'arial',
                              size = 10),
                  textfont = list(textangle = 180))%>%
  
  layout(xaxis = list(title = 'Alignment Position', font='arial'),
         yaxis = list(title = 'Log2 fold difference in mutability', font='arial'))%>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = 400,
      height = 200
    )
  )
p1

# Export the plot with the plotly export (little camera when hovering cursor over plot) for a vector graphic 
# P.S. (add.svg or .pdf manually when saving). The export will have the height and width specified in the "config()" call


# COLORS for PyMOl ################################################################################
# The follow syntax is made for each residue: color 0xffccff, resi 125 and ecLipA_5EXJ_hhpredmodel
# This will be saved in a .pml file and can be easily used to color residues of the pymol model based on
# the calculated fold difference:

# example
color <- "0xffccff"
Res <- "125"
Model_name <- "ecLipA_5EXJ_hhpredmodel"
paste("color ", color, ", ", "resi ", Res, "and ", Model_name, sep = "")

# Loop through the data-frame of residues and their log2 differences, setting a heat-map color for each
# in a new column. The scale for coloring is not linear in a heatmap, so three scales are used depending
# ratio < 0,   0 < ratio < 0.5, and ratio > 0.5


Res_mut$color <- NA
for (R in 1:nrow(Res_mut)){
  # For residues with high log2 ratio values
  if (Res_mut$RatioCount[R] > 0.5){
    # Unfortunately the pymol language uses its own syntax for RBG coloring, so "#" has to be removed and 
    # "0x" has to be added
    Res_mut$color[R] <- paste("0x",substring(rgb(1,1.5-Res_mut$RatioCount[R]*1.5,1.5-Res_mut$RatioCount[R]*1.5), 2),sep="")
  }
  # For residues with low log2 ratio values
  if ((Res_mut$RatioCount[R] < 0.5) & (Res_mut$RatioCount[R] > 0)  ){
    Res_mut$color[R] <- "0xdedede"
  }
  # For residues with negative log2 ratio values
  if (is.na(Res_mut$color[R])){
     #Color everything gray
    Res_mut$color[R] <- "0xdedede"
  }
}
# Make a new column with the PyMol commands for each row/residue
Res_mut$PymolCMD<-paste("color ",Res_mut$color, ", ", "resi ", Res_mut$position, " and ", Model_name, sep = "")
# Save it as a .pml
# When you open this in pymol, it will color the model
write(Res_mut$PymolCMD, "./Structure/pymolcolor.pml")

rgb(1,1.5-(0.5*1.3),1.5-(0.5*1.3))
rgb(1,1.5-(1*1.3),1.5-(1*1.3))

# Plotting a 1 dimensional smoothed heatmap of the overall conservation score for each residue
df <- data.frame(x=Res_mut$position, y=Res_mut$RatioCount)
newdf=data.frame(x=spline(df)[[1]],freq=spline(df)[[2]],y=rep(1,times=3*length(df$x)))

#svg("1dheatmap_all_conservation.svg", height = 1, width = 10)
ggplot(newdf,aes(x=x,y=y,fill=freq))+
  geom_tile()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(fill = "Relative\nresidue\nconservation")
#dev.off()

# Load libraries and data#################################################################################
library(stringr)
library(plyr)
library(ggplot2)
library(seqinr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(ggsci)

# Clear the workspace. set the working dir
rm(list = ls(all.names = TRUE))
dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)

#Reading and cleaning Data   #########################################################################
######################################################################################################

# Read in the fasta file of all LipA sequences from uniprot
All_lipA <- read.fasta(file = "uniprot-name__lipoyl+synthase_.fasta", seqtype = "AA", as.string = TRUE)
# Read the .txt file with ProSiteScan hits with the RSSY motif
#LipA_motifv3<-read.delim("ScanProsite Results for job 210323_V4_3.txt")

# Populated a new data frame with data from motif data
# This takes a while so only do it once
#Motif_data<- data.frame(Name=character(),
                             #Start=integer(), 
                             #End=integer(), 
                            #String=character()) 

#for (n in 1:nrow(LipA_motifv3)){
  #if (substring(LipA_motifv3[n,1], 1,2) == "sp" |substring(LipA_motifv3[n,1], 1,2) == "tr"){
    #df<-data.frame(LipA_motifv3[n,1],LipA_motifv3[n,2],LipA_motifv3[n+1,1],LipA_motifv3[n+6,1])
    #names(df)<-c("Name","Start", "End","String")
    #Motif_data <- rbind(Motif_data, df)
  #}
#}

#write.csv(Motif_data,"Motif_data3.csv")
Motif_data<-read.csv("Motif_data3.csv")


# Make data frame of the fastas of all LipA
Fasta.data<-data.frame(Main=sapply(All_lipA, "[", 1))
Fasta.data$Name <- row.names(Fasta.data)

# Merge the dataframes ###############################
MergedData<-merge(Fasta.data, Motif_data, by = "Name") 
write.csv(MergedData, "./V4 motif split/Motifv4DB_3.csv")

MergedData1<-read.csv("./V4 motif split/Motifv4DB_1.csv")
MergedData2<-read.csv("./V4 motif split/Motifv4DB_2.csv")
MergedData3<-read.csv("./V4 motif split/Motifv4DB_3.csv")
MergedData<-rbind(MergedData1,MergedData2,MergedData3)


# Write all the sequences not covered by a motif, into a fasta file
Fasta.data.nomotif<-Fasta.data[!(Fasta.data$Name %in% MergedData$Name),]
write.fasta(as.list(Fasta.data.nomotif$Main), names=Fasta.data.nomotif$Name, 
            as.string=FALSE, file.out="Not_covered3_v4.fa")
# This was used to make a new DB for Prositescan, for the v4 motif

# Write 100 random no motif sequences for manual inspection of the MSA
#Fasta.data.nomotif.100<-sample_n(Fasta.data.nomotif,100)
#write.fasta(as.list(Fasta.data.nomotif.100$Main), names=Fasta.data.nomotif.100$Name, 
                        #as.string=FALSE, file.out="100RandomNoMotifv4.fa")

# Wokring with and plotting clean data #############################################################
####################################################################################################

# Performed a SAM radical motif scan on the rest of the one not covered
# There seemed to be a lot of sequences that were not actually LipAs, but misanottated as such
# Which will hopefully be caught by this re-annotation
write.fasta(as.list(Fasta.data.nomotif.full$Main), names=Fasta.data.nomotif.full$Name, 
            as.string=FALSE, file.out="Not_covered_byV4motif.fa")


Fasta.data.nomotif$Length<-nchar(Fasta.data.nomotif$Main)
Fasta.data.nomotif.full<-subset(Fasta.data.nomotif,Fasta.data.nomotif$Length > 250)


LipA_motifv3<-read.delim("ScanProsite Results for job 210322 LipA RSSYV4SAM.txt")

# Populated a new data frame with data from motif data
# This takes a while so only do it once
Motif_data<- data.frame(Name=character(),
                        Start=integer(), 
                        End=integer(), 
                        String=character()) 

for (n in 1:nrow(LipA_motifv3)){
  if (substring(LipA_motifv3[n,1], 1,2) == "sp" |substring(LipA_motifv3[n,1], 1,2) == "tr"){
    df<-data.frame(LipA_motifv3[n,1],LipA_motifv3[n,2],LipA_motifv3[n+1,1],LipA_motifv3[n+3,2])
    names(df)<-c("Name","Start", "End","String")
    Motif_data <- rbind(Motif_data, df)
  }
}

### EXPORT THE FINAL FILES
MergedData_final_nomotifV4<-merge(Fasta.data.nomotif, Motif_data, by = "Name") 
write.csv(MergedData_final_nomotifV4, "./V4 motif split/NoMotifV4_withSAMRAD_over250res.csv")
write.fasta(as.list(MergedData_final_nomotifV4$Main), names=MergedData_final_nomotifV4$Name, 
            as.string=FALSE, file.out="NoMotifV4_withSAMRAD_over250res.fa")

### Check for duplicates
MergedData_final_nomotifV4

MergedData_final_nomotifV4_dub <- subset(MergedData_final_nomotifV4,!duplicated(MergedData_final_nomotifV4$Name))

write.fasta(as.list(MergedData_final_nomotifV4_dub$Main), names=MergedData_final_nomotifV4_dub$Name, 
            as.string=FALSE, file.out="NoMotifV4_withSAMRAD_over250res_noDuB.fa")

### Import with phylogenetics
Lip_meta<-read.delim("uniprot-_lipoyl+synthase_ frag")

# Reformat the naming of to be mergeable with the Uniprot output:
for (i in 1:nrow(MergedData_final_nomotifV4)){
  MergedData_final_nomotifV4$Entry.name[i]<-str_split(MergedData_final_nomotifV4$Name[i], "\\|")[[1]][2]
}

colnames(MergedData_final_nomotifV4)<-c("Name", "Sequence","Len", "Start", "End","String", "Entry")
Final_No_V4Motif_SAM<-merge(Lip_meta, MergedData_final_nomotifV4, by = "Entry")


Final_Bac<-subset(Final_No_V4Motif_SAM,Final_No_V4Motif_SAM$Taxonomic.lineage..SUPERKINGDOM.=="Bacteria")
Final_Bac_nofrag<-subset(Final_Bac,Final_Bac$Fragment=="")

ggplot(Final_Bac_nofrag, aes(Taxonomic.lineage..PHYLUM.)) +
  geom_bar(fill = "#0073C2FF")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export the final fasta
write.fasta(as.list(Final_Bac_nofrag$Sequence), names=Final_Bac_nofrag$Entry, 
            as.string=FALSE, file.out="NoMotifV4_withSAMRAD_over250res_noDuB_Bacteria_NoFrag.fa")

# After visual inspection of this alignment of 323 sequences, the following had no conserved R at the conserved position:
NoConsR<-c("A0A650MN45","A0A0E9LLJ2", "A0A1J1LM45", "A0A4R6WXU4", "A0A4V3FVE0", 
  "A0A0K2RDC4", "A0A0K2RDC4", "A0A0S9BYY4", "A0A538D6E7", "F7WZA2", "A0A377UZD9",
  "A0A381H6T7", "A0A2N2UMU0", "A2W9S8", "A0A3D5DAN2", "A0A6M8R671", "A0A3S4GSQ8",
  "A0A5J4KN23", "A0A1J1EUT0", "A0A4Q3J279", "A0A0T5X804", "A0A2E8RUI7", "N9KZX7", 
  "V4Y397", "A0A2E4KC68", "A0A2A2Y8S4", "A0A1Q7TMA7", "A0A2V5YL16", "A0A380E1A1")
Final_Bac_nofrag_noR<-subset(Final_Bac_nofrag, Final_Bac_nofrag$Entry %in% NoConsR)

write.fasta(as.list(Final_Bac_nofrag_noR$Sequence), names=Final_Bac_nofrag_noR$Entry, 
            as.string=FALSE, file.out="Bac_NoR_vis_inspect.fa")

ggplot(Final_Bac_nofrag_noR, aes(Taxonomic.lineage..PHYLUM.)) +
  geom_bar(fill = "#0073C2FF")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



unique(Final_Bac_nofrag_noR$Taxonomic.lineage..PHYLUM.)

noRphyla<-data.frame(Phylum=character(),Total=integer(), NnoR=integer())
noRphyla[1:9,]<-NA

noRphyla$Phylum<-unique(Final_Bac_nofrag_noR$Taxonomic.lineage..PHYLUM.)


for (n in 1:nrow(noRphyla)){
  noRphyla$Total[n]<-nrow(subset(Lip_meta,Lip_meta$Taxonomic.lineage..PHYLUM. == noRphyla$Phylum[n]))
  noRphyla$NnoR[n]<-nrow(subset(Final_Bac_nofrag_noR,Final_Bac_nofrag_noR$Taxonomic.lineage..PHYLUM. == noRphyla$Phylum[n]))
}

noRphyla$Freq<-noRphyla$NnoR/noRphyla$Total

ecname<-"E_coli_LipA"
ecseq<-"MSKPIVMERGVKYRDADKMALIPVKNVATEREALLRKPEWMKIKLPADSTRIQGIKAAMRKNGLHSVCEEASCPNLAECFNHGTATFMILGAICTRRCPFCDVAHGRPVAPDANEPVKLAQTIADMALRYVVITSVDRDDLRDGGAQHFADCITAIREKSPQIKIETLVPDFRGRMDRALDILTATPPDVFNHNLENVPRIYRQVRPGADYNWSLKLLERFKEAHPEIPTKSGLMVGLGETNEEIIEVMRDLRRHGVTMLTLGQYLQPSRHHLPVQRYVSPDEFDEMKAEALAMGFTHAACGPFVRSSYHADLQAKGMEVK"

for (org in Final_Bac_nofrag_noR$Organism){
  orgsub<-subset(Final_Bac_nofrag_noR,Final_Bac_nofrag_noR$Organism==org)
  orgsub[2,5]<-ecname
  orgsub[2,12]<-ecseq
  write.fasta(as.list(orgsub$Sequence), names=c(orgsub$Organism), 
              as.string=FALSE, file.out=paste(org,".fa"))
}



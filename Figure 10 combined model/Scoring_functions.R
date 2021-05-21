# This scripts contains the necessary functions to perform complementation scoring on 
# clustal-omega pairwise alignments of a LipA against E. coli LipA
# David Lennox-Hvenekilde 2021
# Dal@biosyntia.com

# Load libraries
library(seqinr)
library(ggplot2)
library(stringr)
library(tidyr)
library(readr)

# Clear work-space and set wd to the dir where this script is saved
rm(list = ls(all.names = TRUE))
dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)
dir()

# This script is dependent on these input files
# The combined model contains the consveration/complementation data for each residue
# and the epPCR mutability data
Combined_model <- read.csv("Combined_model_rawdatav2.csv")
BLOSUM62<-read.csv("BLOSUM62.csv")
# Clean it up
row.names(BLOSUM62)<-BLOSUM62[,1]
BLOSUM62<-BLOSUM62[,-1]

# read the 15 pairwise alignments of the query gene to ecLipA##############
# ONLY WORK IF ALIGNMENTS HAVE 'align' in the file name and are in a subfolder of the wd()
# You can simply change the pattern from 'align' to another identifier if prefered
# IMPORTANT, the E. coli LipA in the aligments must have the name "E_coli_LipA"
# if not, you can change this identifier to the one you used in the loop below

# List the alignment files (pair-wise alignments)
Align_files <- list.files(pattern = 'align', recursive = TRUE, ignore.case = TRUE)
Align_files[1]

# For each alignment file, this loops through and read them and cleans them and scores them
for ( file in Align_files) {
  temp <- read.alignment(file=file, format = "clustal")
  temp <- as.data.frame(t(as.matrix(temp)))
  temp$EcPOS<-NA
  pos = 1
  for (i in 1:nrow(temp)){
    if ((temp$E_coli_LipA[i] != "\t")&(temp$E_coli_LipA[i]) != "-"){
      temp$EcPOS[i]<-pos
      pos=pos+1
    }
  }
  # Makes a list of which rows to remove (chose with no ecLipA residue)
  row_names_df_to_remove<-c()
  for (i in row.names(temp)){
    if ((temp[i,2]=="\t")|(temp[i,1]=="\t")|(temp[i,2]=="-")|(temp[i,1]=="-")){
      row_names_df_to_remove<-c(row_names_df_to_remove,i)
    }
  }
  # Remove the empty rows where there is no ecLipA residues
  temp<-temp[!(row.names(temp) %in% row_names_df_to_remove),]
  
  # Now we calculate the BLOSUM62 score for each position and store in in the DF
  # For each residue (based in ecLipA) we add the conservation score in a new column
  # In a position where the query sequence or ecLipA sequence is empty then the BLOSUM score = 0
  temp$BLOSUM<-0
  for (pos in 1:nrow(temp)){
    qAA<-str_to_upper(temp[pos,1])
    ECAA<-str_to_upper(temp[pos,2])
    
    if (!(qAA == "-")&!(ECAA == "-") ){
      temp$BLOSUM[pos] <- BLOSUM62[row.names(BLOSUM62)==ECAA,colnames(BLOSUM62)==qAA]
    }
  }
  
  # Now we look through the residues and assign scores
  # The scores are based on the BLOSUM conservation score times the conservation/mutability score:
  temp$Cons_score_RES<-0
  temp$epPCR_score_RES<-0
  for (pos in temp$EcPOS){
    # Calculate and store the Conservation score
    temp[temp$EcPOS==pos,]$Cons_score_RES <-
      Combined_model[Combined_model$EcPOS==pos,]$Compl_cons*
      temp[temp$EcPOS==pos,]$BLOSUM
    
    # Calculate and store the epPCR score score
    temp[temp$EcPOS==pos,]$epPCR_score_RES <-
      Combined_model[Combined_model$EcPOS==pos,]$epPCR_mutability*
      temp[temp$EcPOS==pos,]$BLOSUM
  }
  
  # Change the first column to contain the query gene name
  temp[,1]<-colnames(temp)[1]
  colnames(temp)[1]<-"Query_gene"
  assign(file, temp)
  rm(temp)
}

# Finally we bind all the 15 pairwise alignment files together
# making the final data_fram with which we will 
Final_model_DF<-NA
for (name in Align_files[1:length(Align_files)]) {
  Final_model_DF<-rbind(Final_model_DF,get(name))
}
Final_model_DF<-Final_model_DF[-1,]



# MODELS. See "model_development.R" for details on model generation and refinement 
##################################################################################

OutcomeAlign_SIG = function(Alignment){
  SIGRES_cons<-c(37,78,80,142,174,190,241,273,274,277)
  SIGRES_epPCR<-c(64,66,68,71,73,79,94,96,98,101,108,134,138,139,145,152,171,206,214,218,232,237,240,263,264,265,296,306,308)
  Alignment_sub<-subset(Alignment, Alignment$EcPOS %in% c(SIGRES_cons,SIGRES_epPCR))
  Total=-0.49781+
    0.27414*sum(Alignment_sub$Cons_score_RES^2)+
    0.10395*sum(Alignment_sub$BLOSUM)
  Total
}

OutcomeAlign_SIMPLE = function(Alignment){
  Alignment_sub<-Alignment
  Total=0.275871*sum(Alignment_sub$Cons_score_RES)-0.100733*sum(Alignment_sub$epPCR_score_RES)+0.088921*sum(Alignment_sub$BLOSUM)
  Total
}

# A super simple model, only the sum of epPCR score
OutcomeAlign_epPCR = function(Alignment){
  Alignment_sub<-Alignment
  Total=sum(Alignment_sub$epPCR_score_RES)
  Total
}
# A super simple model, only the sum of conservation score
OutcomeAlign_Cons = function(Alignment){
  Alignment_sub<-Alignment
  Total=sum(Alignment_sub$Cons_score_RES)
  Total
}
# A super simple model, only the sum of conservation score AND epPCR score
OutcomeAlign_epPCR_Cons = function(Alignment){
  Alignment_sub<-Alignment
  Total=sum(Alignment_sub$Cons_score_RES)+sum(Alignment_sub$epPCR_score_RES)
  Total
}
# A super simple model, only the sum of BLOSUM score
OutcomeAlign_BLOSUM = function(Alignment){
  Alignment_sub<-Alignment
  Total=sum(Alignment_sub$BLOSUM)
  Total
}

# Statistical logit models
Scores_SIMPLE<-c()
Scores_SIG<-c()

# Simple additative models
Scores_epPCR<-c()
Scores_epPCR_Cons<-c()
Scores_Cons<-c()
Scores_BLOSUM<-c()

Organisms<-c()


# Look through and score every alignment based on the 6 models
for (A in unique(Final_model_DF$Query_gene)){
  Al_A<-subset(Final_model_DF, Final_model_DF$Query_gene==A)
  
  Scores_SIMPLE<-c(Scores_SIMPLE,OutcomeAlign_SIMPLE(Al_A))
  #Scores_LESSSIMPLE<-c(Scores_LESSSIMPLE,OutcomeAlign_LESSSIMPLE(Al_A))
  #Scores_COMPLEX<-c(Scores_COMPLEX,OutcomeAlign_COMPLEX(Al_A))
  Scores_SIG<-c(Scores_SIG,OutcomeAlign_SIG(Al_A))
  
  Scores_epPCR<-c(Scores_epPCR,OutcomeAlign_epPCR(Al_A))
  Scores_epPCR_Cons<-c(Scores_epPCR_Cons,OutcomeAlign_epPCR_Cons(Al_A))
  Scores_Cons<-c(Scores_Cons,OutcomeAlign_Cons(Al_A))
  Scores_BLOSUM<-c(Scores_BLOSUM,OutcomeAlign_BLOSUM(Al_A))
  
  Organisms<-c(Organisms,A)
  Complementations <- c(Complementations,unique(Al_A$Complementation))
}
# Make a data-frame of the scores
df<- data.frame(Organisms,Scores_SIMPLE,Scores_SIG,
                Scores_epPCR,Scores_epPCR_Cons,Scores_Cons,Scores_BLOSUM)

# You now have a data-frame with the scoring based on the 6 different models
# The Scores_SIG model should give the best prediction

df

##############################################################################################################
# This script is for making a predictive model for LipA complementation based on sequence on a heterologous LipA
# the input is in the form of a pair-wise sequence alignment.
# David Lennox-Hvenekilde 
# 210408

# Load libraries
library(seqinr)
library(ggplot2)
library(stringr)
library(tidyr)
library(readr)


# Clear the workspace. set the working dir. Make sure that this script is in the same folder as the .csv data files
rm(list = ls(all.names = TRUE))
dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)
dir()

# Load the data frames from the two analyses, heterologous LipA and epPCR
epPCR_data <- read.csv("epPCR_Res_mut.csv")
comp_data <- read.csv("Cons_CompvsNoncomp_mtLipA.csv")

# Setting up the alignment for linking data ####################################################
# This alignment is needed as th heterologous complmentation data is based on the mtLipA sequence
# and the epPCR function is based on epLipA sequence. These two need to be connected by the alignment
# Load E coli and M tuberculosis alignment:

alignmentfile <- read.alignment(file="./Alignments/ecLipA_mtLipA_align.clustal_num", format = "clustal")


# Turn the alignment into a data.frame and transpose it, giving a long format
Alignment <- as.data.frame(t(as.matrix(alignmentfile)))

# For some reason there are some empty rows in the alignment data-frame, these are removed:
row_names_df_to_remove<-c()
for (i in row.names(Alignment)){
  if (Alignment[i,2]=="\t"){
    row_names_df_to_remove<-c(row_names_df_to_remove,i)
  }
}
Alignment<-Alignment[!(row.names(Alignment) %in% row_names_df_to_remove),]
Alignment$EcPOS<-1:nrow(Alignment)

# Make a new column and compute the mtLipA residue positions
# This loop takes adds a position (residue number) to any mtLipA row that is not empty
pos = 1
Alignment$MtPOS<-0
for (i in 1:nrow(Alignment)){
  if (Alignment$M_tuberculosis_LipA[i]!="-"){
    Alignment$MtPOS[i]<-pos
    pos = pos + 1
  }
}

# Adding the two model data to the main data alignment file ###########################################
# We now have an alignment data file with the correct positions of the residues in mtLipA and ecLipA
# We can now link the ecLipA (epPCR) and mtLipA (complementation, conservation)

# Add a binary value for the epPCR ration data, same setpoint as used in the original model, just didn't
# add it to the dataframe for some reason
Combined_model <- Alignment

epPCR_data$SIG <- (epPCR_data$RatioCount > -min(epPCR_data$RatioCount))

# For the position in ecLipA, add data on mutability and significance (cutoff) from the epPCR model
Combined_model$epPCR_mutability<-NA
Combined_model$epPCR_mutability_sig<-NA
for (position in Combined_model$EcPOS){
  Combined_model[Combined_model$EcPOS==position,]$epPCR_mutability<-epPCR_data[epPCR_data$position==position,]$RatioCount
  Combined_model[Combined_model$EcPOS==position,]$epPCR_mutability_sig<-epPCR_data[epPCR_data$position==position,]$SIG
}
Combined_model[322,]$epPCR_mutability<-NA
Combined_model[322,]$epPCR_mutability_sig<-NA

# Do the same for mtLipA complementation/conservation and mtLipA
Combined_model$Compl_cons<-NA
Combined_model$Compl_cons_sig<-NA
for (position in Combined_model[Combined_model$MtPOS!=0,]$MtPOS){
  Combined_model[Combined_model$MtPOS==position,]$Compl_cons<-comp_data[comp_data$pos==position,]$DIFF
  Combined_model[Combined_model$MtPOS==position,]$Compl_cons_sig<-comp_data[comp_data$pos==position,]$SIG
}

# Remove the C-terminal of mtLipA, since it does not align to ecLipA
Combined_model[is.na(Combined_model$Compl_cons),]$Compl_cons<-0
Combined_model[is.na(Combined_model$Compl_cons_sig),]$Compl_cons_sig<-0
Combined_model<-subset(Combined_model, !is.na(Combined_model$epPCR_mutability))

# Should remove all negative values for epPCR mutability
Combined_model[Combined_model$epPCR_mutability<0,]$epPCR_mutability<-0
Combined_model[Combined_model$Compl_cons<0,]$Compl_cons<-0
# Export the final file for later use
# Do this only once
#write.csv(Combined_model, "Combined_model_rawdatav2.csv")


# Now we construct the combined model#################################################################
# Clear the workspace and load in the cleaned data
rm(list = ls(all.names = TRUE))
Combined_model <- read.csv("Combined_model_rawdatav2.csv")

BLOSUM62<-read.csv("BLOSUM62.csv")
# Clean it up
row.names(BLOSUM62)<-BLOSUM62[,1]
BLOSUM62<-BLOSUM62[,-1]

# read the 15 pairwise alignments of the query gene to ecLipA##############
# ONLY WORK IF ALIGNMENTS HAVE "align" in the file name and are in a subfolder of the wd()
# List the alignment files (pair-wise alignments)
Align_files <- list.files(pattern = 'align', recursive = TRUE, ignore.case = TRUE)

# For each alignment file, this loops through and read them and cleans them
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

# We need to indicate the outcome, i.e., if the gene complements LipA
unique(Final_model_DF$Query_gene)
complementing<-c("M_tuberculosis_LipA","Pantoea_LipA", "S_liquefaciens_LipA", "P_fluorescens_LipA", 
                 "Synechococcus_sp_LipA", "M_tuberculosis_LipA", "P_torquis_LipA", "P_multisaccharivorax_LipA",
                 "D_teidjei_LipA")
Final_model_DF$Complementation<-0
for (n in 1:nrow(Final_model_DF)){
  if (Final_model_DF$Query_gene[n] %in% complementing){
    Final_model_DF$Complementation[n]<-1
  }
}

#write.csv(Final_model_DF, "Final_model_v2.csv")

# Now we finally have everything we need in one dataframe#######################################
rm(list = ls(all.names = TRUE))
dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)
Combined_model <- read.csv("Combined_model_rawdatav2.csv")
Final_model_DF<-read.csv("Final_model_v2.csv")

# Statistics ###################################################################################
# Here we will make the actual staistical model for complementation as a function of 
# The most significant residues


p <- ggplot(Final_model_DF, aes(x=factor(Query_gene), y=epPCR_score_RES)) + 
  geom_boxplot()
p

t.test(Final_model_DF[Final_model_DF$Complementation==1,]$Cons_score_RES,
       Final_model_DF[Final_model_DF$Complementation==0,]$Cons_score_RES)

hist(Final_model_DF$Cons_score_RES)
hist(Final_model_DF$epPCR_score_RES)

t.test(Final_model_DF[Final_model_DF$Complementation==1,]$epPCR_score_RES,
       Final_model_DF[Final_model_DF$Complementation==0,]$epPCR_score_RES)


# Logit regression #####################################################################
# Make the position a factor, this is needed for the general linear model, otherwise
# it will treat it as a continuous value (i.e. residue number 100 will be 2x more important
# in the model than residue 50 and so forth)


# First we try out a ver simple model that just takes into account the
# Conservation score, and epPCR score (or rather the sum across all residues)
# A general linear model for the logit (logistic regression)
SIGRES_cons<-c(37,78,80,142,174,190,241,273,274,277)
SIGRES_epPCR<-c(64,66,68,71,73,79,94,96,98,101,108,134,138,139,145,152,171,206,214,218,232,237,240,263,264,265,296,306,308)

SIMPLEmodel <- glm(Complementation ~ 
                     Cons_score_RES +
                     epPCR_score_RES +
                     BLOSUM,
                   data = Final_model_DF, family=binomial(link="logit"))
summary(SIMPLEmodel)

# We can drop epPCR from the model
drop1(SIMPLEmodel, test = "Chisq")
# Nothing can be dropped
summary(SIMPLEmodel)

# The model then looks like this
OutcomeAlign_SIMPLE = function(Alignment){
  Alignment_sub<-Alignment
  Total=0.275871*sum(Alignment_sub$Cons_score_RES)-0.100733*sum(Alignment_sub$epPCR_score_RES)+0.088921*sum(Alignment_sub$BLOSUM)
  Total
}

# A less simple model that takes into acount the interactions between conservation score
# and epPCR mutability
# We do not inlude interactions between BLOSUM and epPCR/cons as they are already products
LESSSIMPLEmodel <- glm(Complementation ~ 
                         Cons_score_RES +
                         I(Cons_score_RES^2)+
                         epPCR_score_RES +
                         I(epPCR_score_RES^2)+
                         I(Cons_score_RES*epPCR_score_RES)+
                         BLOSUM+
                         I(BLOSUM^2),
                       data = Final_model_DF, family=binomial(link="logit"))
summary(LESSSIMPLEmodel)

drop1(LESSSIMPLEmodel, test = "Chisq")
LESSSIMPLEmodel <- update(LESSSIMPLEmodel,~.-I(Cons_score_RES^2))
drop1(LESSSIMPLEmodel, test = "Chisq")
LESSSIMPLEmodel <- update(LESSSIMPLEmodel,~.-I(Cons_score_RES * epPCR_score_RES))
drop1(LESSSIMPLEmodel, test = "Chisq")
LESSSIMPLEmodel <- update(LESSSIMPLEmodel,~.-I(epPCR_score_RES^2))
drop1(LESSSIMPLEmodel, test = "Chisq")
LESSSIMPLEmodel <- update(LESSSIMPLEmodel,~.-I(BLOSUM^2))
summary(LESSSIMPLEmodel)
# This less simple model is reduced to the simple model

# Now we make a general linear model for the logit (logistic regression)
# this is a statistical way of modeling a binomial (complementation) outcome with one or 
# more explanatory variables (the conservation and epPCR mutability score at each residue in alignment
# rooted in the ecLipA sequence residue number)


# TAKE A LONG TIME AND CPU TO RUN THE MODEL, DO NOT RUN IT UNLESS NECESSARY
# Load the already run model instead (below)
Final_model_DF_SIG<-subset(Final_model_DF, Final_model_DF$EcPOS %in% c(SIGRES_cons,SIGRES_epPCR))


mylogitSIG <- glm(Complementation ~ 
                 Cons_score_RES +
                 I(Cons_score_RES^2)+
                 epPCR_score_RES +
                 I(epPCR_score_RES^2)+
                 I(Cons_score_RES*epPCR_score_RES)+
                 BLOSUM+
                 I(BLOSUM^2),
                 data = Final_model_DF_SIG, family=binomial(link="logit"),maxit = 100)
summary(mylogitSIG)

drop1(mylogitSIG, test = "Chisq")
mylogitSIG <- update(mylogitSIG,~.-Cons_score_RES)
drop1(mylogitSIG, test = "Chisq")
mylogitSIG <- update(mylogitSIG,~.-epPCR_score_RES)
drop1(mylogitSIG, test = "Chisq")
mylogitSIG <- update(mylogitSIG,~.-I(epPCR_score_RES^2))
drop1(mylogitSIG, test = "Chisq")
mylogitSIG <- update(mylogitSIG,~.-I(Cons_score_RES * epPCR_score_RES))
drop1(mylogitSIG, test = "Chisq")
mylogitSIG <- update(mylogitSIG,~.-I(BLOSUM^2))
drop1(mylogitSIG, test = "Chisq")
summary(mylogitSIG)


OutcomeAlign_SIG = function(Alignment){
  Alignment_sub<-subset(Alignment, Alignment$EcPOS %in% c(SIGRES_cons,SIGRES_epPCR))
  Total=-0.49781+
    0.27414*sum(Alignment_sub$Cons_score_RES^2)+
    0.10395*sum(Alignment_sub$BLOSUM)
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

# Lets test out these models on all the alignments we have:##################################

# Statistical models
Scores_SIMPLE<-c()
Scores_SIG<-c()

# Simple additative models
Scores_epPCR<-c()
Scores_epPCR_Cons<-c()
Scores_Cons<-c()
Scores_BLOSUM<-c()

Organisms<-c()
Complementations<-c()
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
df<- data.frame(Organisms,Complementations,Scores_SIMPLE,Scores_SIG,
                Scores_epPCR,Scores_epPCR_Cons,Scores_Cons,Scores_BLOSUM)
df
write.csv(df,"Final_complementation_scores_logistical_regressionv2.csv")
# It seems to predict pretty well, except one outlier, which is Synechococcus LipA for some reason

# Look at the 8 prediction models#############################################################
rm(list = ls(all.names = TRUE))
models<-read.csv("Final_complementation_scores_logistical_regressionv2.csv")

# long format for plotting
models_long <- models %>%  gather(Model,Score,Scores_SIMPLE:Scores_BLOSUM )


p <- ggplot(models_long, aes(x=Model, y=Score, fill=factor(Complementations))) + 
  geom_boxplot()+
  facet_wrap(~Model, scale="free")
p
p+ scale_fill_manual(values=c("#E64B35", "#00A087"))


PredictionModel<-c()
pvalue<-c()

for (model_type in unique(models_long$Model)){
  PredictionModel<-c(PredictionModel,model_type)
  pvalue<-c(pvalue,as.numeric(t.test(subset(models_long,(models_long$Complementations==0) & (models_long$Model==model_type))$Score,
                                     subset(models_long,(models_long$Complementations==1) & (models_long$Model==model_type))$Score)[3]))
}
AllPvalues <-data.frame(PredictionModel,pvalue)
write.csv(AllPvalues,"Pvalues_for_models.csv")





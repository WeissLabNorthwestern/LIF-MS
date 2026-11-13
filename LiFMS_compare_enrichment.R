#!/usr/bin/env Rscript
#setwd('/home/weisslab/Documents/Ben/Notebook/260611-mess_with_lifms_nm2/all_oxidations/nkd')
# setwd('/home/weisslab/Documents/Ben/Notebook/260721-mess_with_lifms_hmob1b/')
# Score crosslinked residues against total score for each residue; IE,
# how do the scores of the crosslinked and uncrosslinked residues compare?
# Benjamin Parker, PhD
# 6.15.26

# Usage: first command is filename. Second is the prefix (FDR or noFDR) depending on which label extractor
# version is used-one cuts off the dataset at the highest ranking DECOY peptide, the other doesn't.
# Third is the minimum score value.

ExperimentID <- commandArgs(trailingOnly = TRUE)
#DEBUG
# ExperimentID <-'mob1b-eddei-cterm.txt'
# data_prefix<-ExperimentID
# is_FDR<-'FDR'
data_prefix<-ExperimentID[1]
# Which FDR output suffix for filename?
is_FDR<-ExperimentID[2]
# Minimum score to export? 
min_score<-ExperimentID[3]

# Read in files.
tagged<-read.csv(paste(data_prefix,'-labels.csv-',is_FDR,sep=''))
total<-read.csv(paste(data_prefix,'-total_score.csv',sep=''))

names(total)<-c("butanolized_residue","butanol_location","total_score")

# Merge together, matching the crosslinked residues in the protein sequence with the uncrosslinked ones.
all_residues_scored<-merge(total, tagged,all.y=TRUE)
all_residues_scored[is.na(all_residues_scored)] <- 0
#all_residues_scored$tag_score_fraction<-(all_residues_scored$score/all_residues_scored$total_score*20144623751753) # shift decimal point by 14 digit prime number
all_residues_scored$tag_score_fraction<-(all_residues_scored$score/log10(all_residues_scored$total_score)) # shift decimal point by 14 digit prime number

all_residues_scored[is.infinite(all_residues_scored$tag_score_fraction),5] <-100003# 6 digit prime number 102221243011	# Replace infinite values with a 12-digit prime number

# Filter out  scores below min_score
all_residues_scored <-all_residues_scored[all_residues_scored$tag_score_fraction>min_score,]

# Final, exported data frame
scores_toexport <- all_residues_scored[,c(1,2,5)]
names(scores_toexport) <-c("butanolized_residue","butanol_location","score")
write.table(scores_toexport[order(scores_toexport$score,decreasing=TRUE),], paste(ExperimentID[1],'-normalized_score',sep=''), col.names=TRUE, sep=",", row.names=FALSE)

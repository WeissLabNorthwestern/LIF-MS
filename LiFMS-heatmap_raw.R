#!/usr/bin/env Rscript
library(bio3d)
suppressMessages(library(dplyr))
# Convert LIFMS input and input pdb into a pdb file with crosslink
# Scores replacing the B-factor values.
# Benjamin Parker, Ph.D.
# June 1st, 2026

# Modified August 5th, 2026o to remove ranking scoring system.

 ExperimentID <- commandArgs(trailingOnly = TRUE)
 score_file_name <-ExperimentID[1]
 pdb_file_name <- ExperimentID[2]
# 
# # #DEBUG
#setwd('/home/weisslab/Documents/Ben/Notebook/260805-lifms_ndr_final_pipeline/')
# score_file_name<-'nkd_sample1-2-bidirectionalized_order'
# score_file_name<-'nkd_sample3-4-bidirectionalized_order'
# score_file_name<-'nkd_sample5-6-bidirectionalized_order'
#score_file_name<-'nkd_sample7-8-bidirectionalized_order'

# 
#pdb_file_name<-'6bxi-cleaned.pdb'
# # 


input_PDB <-read.pdb(pdb_file_name)
input_scores_raw<-read.csv(score_file_name)


# input_scores_raw$order_normalized<-input_scores_raw$order/max(input_scores_raw$order)
input_scores<-input_scores_raw |> select(butanolized_residue,butanol_location,order_normalized)
names(input_scores) <-c("butanolized_residue","butanol_location","score")

# Extract atoms from data frame. We need to merge the crosslinks with the atom values,
# Then re-input the crosslinks and scores back into the pdb file.
input_atoms<-as.data.frame(input_PDB$atom$resno)
names(input_atoms)<-('butanol_location')
merged_scores_raw <-merge(input_atoms,input_scores,all.x=TRUE)

merged_scores<-merged_scores_raw
# Remove butanol locations that are not actual residues (IE, from decoys.)
# merged_scores <-merged_scores_raw[which(!merged_scores_raw$butanol_location==0),]


# Compress the scores-they have a massive dynamic range that is hard to visualize.
merged_scores[is.na(merged_scores)] <-0



# Replace 0 with 1, since will take

# Form output pdb, replacing current B-factors with the LIFMS scores.
output_pdb <-input_PDB
output_pdb$atom$b <-merged_scores$score

write.pdb(output_pdb,paste(score_file_name,'-',pdb_file_name,'-heatmap','.pdb',sep=''))

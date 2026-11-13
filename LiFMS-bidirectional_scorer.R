#!/usr/bin/env Rscript
library(bio3d)
suppressMessages(library(dplyr))
# Score two-ended LiFMS datasets. Use output from LiFMS_compare_enrichment.R
# August 5th, 2026
# Benjamin Parker, PhD
setwd('/home/steven/Documents/Summary_LIFMS_work_Ndr/LSSSP/comet/')

# Load in data frames.
nterm_data <- read.csv('nkd-sample1.txt-normalized_score');cterm_data <- read.csv('nkd-sample2.txt-normalized_score');export_filename <-'nkd_sample1-2'
# nterm_data <- read.csv('nkd-sample3.txt-normalized_score'); cterm_data <- read.csv('nkd-sample4.txt-normalized_score');export_filename <-'nkd_sample3-4'
# nterm_data <- read.csv('nkd-sample5.txt-normalized_score');cterm_data <- read.csv('nkd-sample6.txt-normalized_score');export_filename <-'nkd_sample5-6'
# nterm_data <- read.csv('nkd-sample7.txt-normalized_score');cterm_data <- read.csv('nkd-sample8.txt-normalized_score');export_filename <-'nkd_sample7-8'
# nterm_data <- read.csv('nkd-sample9.txt-normalized_score');cterm_data <- read.csv('nkd-sample10.txt-normalized_score');export_filename <-'nkd_sample9-10'
# nterm_data <- read.csv('nkd-sample11.txt-normalized_score');cterm_data <- read.csv('nkd-sample12.txt-normalized_score');export_filename <-'nkd_sample11-12'

# 
# 

# 

# Off target crosslinks to remove.
garbage_sites<-c(seq(82,86),seq(122,127),seq(259,264))

# Order crosslinked residues.
nterm_data <-nterm_data[order(nterm_data$score, decreasing=FALSE),]
nterm_data$order <- seq(1:nrow(nterm_data))
nterm_data$order_normalized <-nterm_data$order/max(nterm_data$order)
nterm_data$score_normalized <-nterm_data$score/max(nterm_data$score)

nterm_data_final<-nterm_data |> select(butanolized_residue,butanol_location,score,order_normalized)
# nterm_data$score<-log10(nterm_data$score)

cterm_data <-cterm_data[order(cterm_data$score, decreasing=FALSE),]
cterm_data$order <- seq(1:nrow(cterm_data))
cterm_data$order_normalized <-cterm_data$order/max(cterm_data$order)
cterm_data$score_normalized <-cterm_data$score/max(cterm_data$score)
cterm_data_final<-cterm_data |> select(butanolized_residue,butanol_location,score,order_normalized)



# Merge, remove garbage sites, then subtract score of Cterm from score of Nterm.
merged_data_raw<-merge(nterm_data_final,cterm_data_final,by=c('butanolized_residue','butanol_location'),all=TRUE)
merged_data <-merged_data_raw[!merged_data_raw$butanol_location %in% garbage_sites,]
# merged_data<-merged_data_raw
names(merged_data)<-c('butanolized_residue','butanol_location','score_nterm','order_nterm','score_cterm','order_cterm')
merged_data[is.na(merged_data)]<-0

merged_data$score_nterm_normalized <-merged_data$score_nterm/max(merged_data$score_nterm)
merged_data$score_cterm_normalized <-merged_data$score_cterm/max(merged_data$score_cterm)

merged_data$n_minus_c_diff <-merged_data$score_nterm_normalized^(1/5)-merged_data$score_cterm_normalized^(1/5)




# Split data, then "re-order" by rank. In this case, the order is by difference between
# the two peptides, not the actual score.
# merged_data_nterm<-merged_data[merged_data$n_minus_c_diff>0,]
# merged_data_nterm <-merged_data_nterm[order(merged_data_nterm$n_minus_c_diff, decreasing=FALSE),]
# merged_data_nterm$order <- seq(1:nrow(merged_data_nterm))
# merged_data_nterm$order_normalized <-merged_data_nterm$order/max(merged_data_nterm$order)
# # Cterm
# merged_data_cterm<-merged_data[merged_data$n_minus_c_diff<0,]
# merged_data_cterm <-merged_data_cterm[order(merged_data_cterm$n_minus_c_diff, decreasing=TRUE),]
# merged_data_cterm$order <- seq(1:nrow(merged_data_cterm))
# # Negative, since negative is C-term.
# merged_data_cterm$order_normalized <--merged_data_cterm$order/max(merged_data_cterm$order)
# 
# # Concatenate the data frames.
# merged_data_toexport<-rbind(merged_data_nterm
#                             |> select(butanolized_residue,butanol_location,order_normalized),
#                             merged_data_cterm
#                             |> select(butanolized_residue,butanol_location,order_normalized)
# )


# Try normalizing the positive and negative differences between scores to 1.
# #
# merged_data_nterm<-merged_data[merged_data$n_minus_c_diff>0,]
# merged_data_nterm <-merged_data_nterm[order(merged_data_nterm$n_minus_c_diff, decreasing=FALSE),]
# merged_data_nterm$normalized_score <-merged_data_nterm$n_minus_c_diff/max(merged_data_nterm$n_minus_c_diff)
# # Cterm
# merged_data_cterm<-merged_data[merged_data$n_minus_c_diff<0,]
# merged_data_cterm <-merged_data_cterm[order(merged_data_cterm$n_minus_c_diff, decreasing=TRUE),]
# merged_data_cterm$normalized_score <--merged_data_cterm$n_minus_c_diff/min(merged_data_cterm$n_minus_c_diff)
# 
# 
# # Concatenate the data frames.
# merged_data_toexport<-rbind(merged_data_nterm
#                             |> select(butanolized_residue,butanol_location,normalized_score),
#                             merged_data_cterm
#                             |> select(butanolized_residue,butanol_location,normalized_score)
# )



# merged_data_toexport <- merged_data |> select(butanolized_residue,butanol_location,normalized_score)
# 
merged_data_toexport <- merged_data |> select(butanolized_residue,butanol_location,n_minus_c_diff)

merged_data_toexport$n_minus_c_diff <-merged_data_toexport$n_minus_c_diff+2

names(merged_data_toexport)<-c('butanolized_residue','butanol_location','order_normalized')



write.table(merged_data_toexport[order(merged_data_toexport$order_normalized,decreasing=TRUE),], paste(export_filename,'-bidirectionalized_order',sep=''), col.names=TRUE, sep=",", row.names=FALSE)

#!/usr/bin/env Rscript

library(bio3d)
library(cluster)
library(stringr)

suppressMessages(library(dplyr))
# LIFMS Cluster junkie
# Benjamin Parker, Ph.D.
# 2.18.2026
#
# Cluster crosslinked resides from the label extractor distance_lifms_script
# on a pdb file. Then, hierarchically cluster the residues based on distance.
#
# Each cluster is a binding site and is scored accordingly by multiplying
# the scores of each residue in the cluster togetheer.





# Name of input file
ExperimentID <- commandArgs(trailingOnly = TRUE)
# Name of output file.
OutputFileBOH <- paste(ExperimentID[1],"-clusters",sep ='')

# Read input file. Remove NAs.
residues_raw<-read.csv(ExperimentID[1])
residues<-residues_raw[!is.na(residues_raw$butanolized_residue),]

# Read structure file.
ref_structure <- read.pdb(ExperimentID[2])

# Create a table of the locations of each atom of the pdb found in the
# crosslinked ("butanol") dataset.
location_table <-ref_structure$atom$resno %in% residues$butanol_location

# Filter atoms from the pdb by whether the are found in the crosslinked dataset.
tagged_residues<-ref_structure$atom[location_table,]

### Keep the C alphas, since these will be used to measure distances.
# Make a data frame of the XYZ coordinate of each C alpha from each crosslinked residue.
tagged_residues$calpha<-ref_structure$calpha[location_table]
tagged_ca<-tagged_residues[tagged_residues$calpha==TRUE,]
# tagged_ca_coords<-as.data.frame(NA)
tagged_ca_coords<-as.data.frame(tagged_ca$resno)
tagged_ca_coords$x<-tagged_ca$x
tagged_ca_coords$y<-tagged_ca$y
tagged_ca_coords$z<-tagged_ca$z
names(tagged_ca_coords)<-c('resno','x','y','z')

### Put the score from the crosslinked residue in the C alpha data frame
tagged_ca_coords$score<-residues[residues$butanol_location %in% tagged_ca_coords$resno,]$score
names(tagged_ca_coords)<-c('resno','x','y','z','score')


### Create a distance matrix of NxN rows, with 1..N being the number of crosslinked residue.
# Iterate from each residue and measure the magnitude of the distance between
# the C alpha from each residue to each other residue.
# The resulting matrix is the distance from each crosslinked residues' C alpha to each
# other residues' C alpha.
rmsd_table_raw<-matrix(nrow=nrow(tagged_ca_coords),ncol=nrow(tagged_ca_coords))
for(i in 1:nrow(tagged_ca_coords)){
  for(j in 1:nrow(tagged_ca_coords)){
    calpha_i<-as.numeric(tagged_ca_coords[i,2:4])
    calpha_j<-as.numeric(tagged_ca_coords[j,2:4])
    subtracted_coords<-calpha_i-calpha_j
    vector_mag<-sqrt(subtracted_coords[1]^2+subtracted_coords[2]^2+subtracted_coords[3]^2)
    rmsd_table_raw[i,j]<-vector_mag

}
}

### Hierarchically cluster the residues.
# The number of clusters is determined by
# the rounded standard deviation of the distance
# matrix divided by 2. This gives the best results in the
# control LIFMS datasets for some reason.
d <- dist(rmsd_table_raw, method = "manhattan")
hc <- agnes(d, method = "complete" )
sub_hc <- cutree(hc, k = round(sd(rmsd_table_raw)/exp(1)))



### Create data table from distance matrix.
#  The distance matrix for each residue is summed for each residue.
# This is effectively its "global distance" from the others.
# This isn't really used in this version of the script.
rmsd_table_raw[is.infinite(rmsd_table_raw)]<-0
rmsd_table_df <-as.data.frame(tagged_ca_coords$resno)
rmsd_table_df$av_rmsd<-apply(rmsd_table_raw,1,sum)
rmsd_table_df$cluster_id<-sub_hc
names(rmsd_table_df)<-c('butanol_location','av_rmsd','cluster_id')
rmsd_table_df_final <-merge(rmsd_table_df,residues_raw)
# rmsd_table_df_final$geometry_score <-rmsd_table_df_final$score*rmsd_table_df$av_rmsd
# 
# write.table(clusters, OutputFileBOH, col.names=TRUE, sep=",", row.names=FALSE,quote=FALSE)

# Multiply
cluster_summary<-rmsd_table_df_final %>% group_by(cluster_id) %>% summarize_at(vars(score),prod)
names(cluster_summary)<-c('cluster_id','cluster_score')
rmsd_table_df_final_clustered<-merge(cluster_summary,rmsd_table_df_final)
rmsd_table_df_final_clustered$cluster_score<-log10(rmsd_table_df_final_clustered$cluster_score)
write.table(rmsd_table_df_final_clustered,OutputFileBOH , col.names=TRUE, sep=",", 
            row.names=FALSE,quote=FALSE)

### Output selectable string for pymol ###
resi_string<-rmsd_table_df_final_clustered[order(rmsd_table_df_final_clustered$cluster_score,decreasing=TRUE),]

# Print the top cluster's residues for pymol.
resi_string_final <-paste(resi_string$butanol_location[ which(resi_string$cluster_id==resi_string$cluster_id[1])],collapse="+")
cat(paste(ExperimentID[1],",",resi_string_final,"\n",sep=""),file="Binding_sites.txt",append = TRUE)

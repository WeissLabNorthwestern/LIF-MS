#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(plyr)

# set filename
input_header = commandArgs(trailingOnly = TRUE)
#input_header = 'MKK4-PepI_onto_JNK1.txt-butanols.csv'
# Input library
boh_input <-read.csv(input_header,skip=1,header=FALSE)
# Currently set for length of human stk38 (ndr1).
protein_length <- 400

# Create data frame of total residue positions for protein.
ref <- data.frame(c(1:protein_length))

# Rename headers.
names(boh_input) <-c('Residue','Number','NLEV')
names(ref) <-c('Number')


# Remove NAs
boh_input_noNA <- boh_input[is.na(boh_input$Residue) == FALSE,]

# Sort by residue number.
boh_input_sorted <- boh_input_noNA[order(boh_input_noNA$Number),]

# Remove NAs, if present.
boh_input_sorted$NLEV[is.na(boh_input_sorted$butanolized_residue) == TRUE] <- 0

# Convert residue number to character so it doesn't treat it as a number.
boh_input_sorted$Number <- as.character(boh_input_sorted$Number)
#boh_input_sorted$Threshold <- as.character(ifelse(boh_input_sorted$NLEV >= 5,1,0))

# Merge data frames with reference.
boh_allres <- merge(ref,boh_input_sorted,by.x = 'Number',by.y = 'Number',all.x = TRUE,all.y = TRUE)

# Fill missing spotes in merged data frame.
boh_allres[is.na(boh_allres)] <- 0

# Graph the results.
ggplot(boh_allres,aes(as.numeric(Number),log10(NLEV))) + 
  geom_bar(stat="identity",aes(fill=NLEV)) +
  scale_fill_gradient(low='darkgrey', high='green') +
scale_colour_gradient(trans = "log10")+
theme(legend.position = 'none',
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      text = element_text(size=12))+ 
   scale_x_discrete(breaks=seq(0,protein_length,100)) +
  xlim(0,protein_length) +
  theme(panel.background = element_rect(fill="aliceblue"))
  
  
ggsave(paste(input_header,"-graph.png",sep=""),width=12,height = 2.5,units = "cm")
# width was 24 before 8.13.19

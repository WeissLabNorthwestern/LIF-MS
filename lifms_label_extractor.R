#!/usr/bin/env Rscript

# Extract labels from a tsv of tryptic peptides; E.G. from COMET. 
suppressMessages(library(dplyr))
# dplyr version 0.7.4 is used.
library(stringr)
# stringr version 1.2.0 is used. 

# Benjamin Parker, PhD
# UPDATES 4.15.24
# Can remove QN as well as IAA-alkylated cys.

# ##### USAGE###
# lifms_label_extractor.R [txt file from comet] [mass tag to seach for as listed in comet] [location of reference fasta]
# For output.txt, butylpyrrolidine mass tag (125.1204Da) and the STK38 kinase, 
# an example command might be:
#  lifms_label_extractor.R output.txt 125.1204 stk38.fasta

# Outputs a list of crosslink tag locations based on the residues from the input fasta. Make
# sure that COMET is using the same fasta as this script when generating the txt of tryptic peptides.
# FDR is calculated by first location of decoy in list of peptides; see the noFDR script if you don't want this.


### REVISED FOR NEW COMET 2018 VERSION ###
# Comet no longer uses rev_sp for decoy, instead using the DECOY prefix.

### Butanol extractor script. 
# Script is designed to take an input from comet MS/MS and determine the
# location in a given protein sequence of each butanol modification. 
# Output .csv contains residues modified by butanol as found in the
# fasta sequence. Ensure the same sequence is used in this file as
# was found in the fasta reference sequence file used by comet. 




# File from Comet without the .txt ending.
# The default settings assumes the file is produced
# by comet with no other modifications. 
#setwd("/home/weisslab/Documents/Ben/Notebook/240318-Compare_NM2_LIFMS_data_rollnotebook17page10/2024_data")
#ExperimentID <- "CGLSS_pyrl_lysis_EXP.txt"
ExperimentID <- commandArgs(trailingOnly = TRUE)

# What is the molecular weight of the tag (IE, butanol, butylamine) you are looking for?
 #tagmass<-71.0735 #Butylamine
 #tagmass <-"72.0575" #Butanol
 #tagmass <-98.0724
#tagmass <- "125.1204" #Butylpyrrolidine-2H
#tagmass<-"128.0956" # iodoacetamide alylkated butylamine
#tagmass <-71.0503 #Butanol with no proton due to CID of crosslinlker
#tagmass <-"125.1215"

# 4.15.24
# Input mass tag mass as command line.

tagmass <-as.character(ExperimentID[2])

# 5.2.24
# Input protein sequence as fasta file
protein_sequence_file_location <-ExperimentID[3]
#protein_sequence_file_location <-"/home/weisslab/Documents/scriptcore/human_proteome/p38a.fasta"

# File name variables. Change as necessary.
# This adds the file name to an extension string. Comet outputs
# .txt, but we usually output .csv or .tsv.
# This ensures your ExperimentID is the same in all casse.
# OutputFileBOH is a .csv containing locations of butanol-modified residues in the 
# protein amino acid sequence.
# OutputFilePeps is a .csv containing butanol-modified peptides. 
# Other modifications (methionine oxidations, deamidations of N and Q, etc)
# have been removed from this data, but butanol modifications remain.  
#InputFile <- paste(ExperimentID,'.txt',sep='')
#setwd('/home/weisslab/Documents/Ben/Notebook/240318-Compare_NM2_LIFMS_data_rollnotebook17page10/2019_data')
#ExperimentID <- "ndr_lifms01-mystery128"
#InputFile <- paste(ExperimentID[1],'.txt',sep='')
InputFile <- ExperimentID[1] # Removed paste 5.15.24
OutputFileBOH <- paste(ExperimentID[1],"-labels.csv",sep ='')
OutputFilePeps <- paste(ExperimentID[1],"-peptides.csv",sep ='')


# Read in .txt PSM list from comet. This reads in a csv with a tab separation instead of a comma.
# To change to comma separation, remove sep ="\t", which specifies tab separation.
# Comet usually adds a header (comet version, protein name, etc) which is only one column; 
# this conflicts with R, so skip=1 for files freshly made from comet to ensure the read function skips it. 
# If your file doesn't have a header like this, set skip=0.

Sample_raw <- read.csv(InputFile,sep ="\t",row.names = NULL,skip=1)

# Confidence cutoff (set to 0). 
confidence_cutoff <- 0 # 4.21.24; was 1

# E-value cutoff
e_val_cutoff=1

# The swissprot FASTA protein amino acid sequence. See examples-it needs to be the swissprot format.
# Note that the sequence must be in QUOTES (""), with the first set of quotes on the same
# line as the FASTA header (see example files).
# The script reads in the header, looking for the first newline, and assumes
# everything on the subsequence lines is the protein amino acid sequence.

# The header also is filtered to extract the fasta entry_name (for example, MK08_HUMAN), so ensure
# the sequence as a properly formatted swissprot format FASTA header.
#  

protein_sequence_INPUT <- readLines(protein_sequence_file_location)
#print(protein_sequence_INPUT)
# protein_sequence_INPUT <- 
# ">sp|Q15208|STK38_HUMAN Serine/threonine-protein kinase 38 OS=Homo sapiens OX=9606 GN=STK38 PE=1 SV=1
# MAMTGSTPCSSMSNHTKERVTMTKVTLENFYSNLIAQHEEREMRQKKLEKVMEEEGLKDE
# EKRLRRSAHARKETEFLRLKRTRLGLEDFESLKVIGRGAFGEVRLVQKKDTGHVYAMKIL
# RKADMLEKEQVGHIRAERDILVEADSLWVVKMFYSFQDKLNLYLIMEFLPGGDMMTLLMK
# KDTLTEEETQFYIAETVLAIDSIHQLGFIHRDIKPDNLLLDSKGHVKLSDFGLCTGLKKA
# HRTEFYRNLNHSLPSDFTFQNMNSKRKAETWKRNRRQLAFSTVGTPDYIAPEVFMQTGYN
# KLCDWWSLGVIMYEMLIGYPPFCSETPQETYKKVMNWKETLTFPPEVPISEKAKDLILRF
# CCEWEHRIGAPGVEEIKSNSFFEGVDWEHIRERPAAISIEIKSIDDTSNFDEFPESDILK
# PTVATSNHPETDYKNKDWVFINYTYKRFEGLTARGAIPSYMKAAK
# "

#protein_sequence_INPUT <- 
#">sp|Q16539|MK14_HUMAN Mitogen-activated protein kinase 14 OS=Homo sapiens OX=9606 GN=MAPK14 PE=1 SV=3
#MSQERPTFYRQELNKTIWEVPERYQNLSPVGSGAYGSVCAAFDTKTGLRVAVKKLSRPFQ
#SIIHAKRTYRELRLLKHMKHENVIGLLDVFTPARSLEEFNDVYLVTHLMGADLNNIVKCQ
#KLTDDHVQFLIYQILRGLKYIHSADIIHRDLKPSNLAVNEDCELKILDFGLARHTDDEMT
#GYVATRWYRAPEIMLNWMHYNQTVDIWSVGCIMAELLTGRTLFPGTDHIDQLKLILRLVG
#TPGAELLKKISSESARNYIQSLTQMPKMNFANVFIGANPLAVDLLEKMLVLDSDKRITAA
#QALAHAYFAQYHDPDDEPVADPYDQSFESRDLLIDEWKSLTYDEVISFVPPPLDQEEMES
#"

#protein_sequence_INPUT <- 
#">sp|Q70IA6|MOB2_HUMAN MOB kinase activator 2 OS=Homo sapiens OX=9606 GN=MOB2 PE=1 SV=1
#MDWLMGKSKAKPNGKKPAAEERKAYLEPEHTKARITDFQFKELVVLPREIDLNEWLASNT
#TTFFHHINLQYSTISEFCTGETCQTMAVCNTQYYWYDERGKKVKCTAPQYVDFVMSSVQK
#LVTDEDVFPTKYGREFPSSFESLVRKICRHLFHVLAHIYWAHFKETLALELHGHLNTLYV
#HFILFAREFNLLDPKETAIMDDLTEVLCSGAGGVHSGGSGDGAGSGGPGAQNHVKER
#"


####### BEGIN SCRIPT ####### 

# Command to extract the protein entry_name from the FASTA file header above.
# (For example, JNK1 is MK08_HUMAN, p38 is MK14_HUMAN, etc).
# 5.2.24-needs to be first line.
Protein_code <- str_extract(
  str_extract(
    protein_sequence_INPUT[1],"(?<=\\|)(.*?)(?=\\s)"
    )
  ,"(?<=\\|)(.*)"
  )

# Everything before the first newline character and the
# ">" (with the FASTA header) is discarded by these str_sub() commands.
# The remaining lines are the amino acid sequence of the protein interspersed with newline \n 
# characters, which are removed by str_replace_all().
# Final protein_sequence variable is the pure protein sequence with only amino acids. 

protein_sequence <-str_sub(
  protein_sequence_INPUT[2],regexpr(">", protein_sequence_INPUT)
  )
# 5.2.24
# All non-fasta lines (protein sequence)
protein_sequence <-  protein_sequence_INPUT[2:length(protein_sequence_INPUT)]
# 
# protein_sequence <- str_replace_all(
#   str_sub(
#     protein_sequence,
#           regexpr("\\n", protein_sequence)
#           ), 
#   "[^a-zA-Z]", ""
#   )
#print(protein_sequence)
# Rename output .csv column names, as R tries to change them.
names(Sample_raw) <- c("scan",	"num",	"charge",	"exp_neutral_mass",	"calc_neutral_mass",	
"e-value",	"xcorr",	"delta_cn",	"sp_score",	"ions_matched",	"ions_total",	"plain_peptide",	
"modified_peptide",	"prev_aa",	"next_aa",	"protein",	"protein_count",	"modifications","blank") 

# FDR CUTOFF
# Sort in order of increasing e-value
sortedraw <-Sample_raw[order(Sample_raw$`e-value`),]

# Vector containing locations in sorted data of reverse peptides.
# 4.15.24; if no DECOY found, set limit to 1000.
#if ((isTRUE(grep('DECOY',sortedraw$protein)))){
#decoy_locations <- grep('DECOY',sortedraw$protein)
#}else{
#decoy_locations <-c('1000')
#}
# 4.15.24; reverting to old command due to above if statement not working
# Don't feel like debuggin it.
decoy_locations <- grep('DECOY',sortedraw$protein)
#decoy_locations <- grep(999,sortedraw$3-value) # Ben test 5.24.24

# Location of the first decoy (this establishes the FDR.)
fdrpoint=decoy_locations[1]

# Raw data from the start to the first decoy (IE, the highest e-valued decoy.)
fdrset <- sortedraw[1:fdrpoint,]


# Filter out PSMs containing butanols.
# Butylamine
#Marked_peptides <- dplyr::filter(fdrset,grepl("71.0735",modified_peptide))
# Butyl azide
Marked_peptides <- fdrset[grep(tagmass,fdrset$modified_peptide),]
 # dplyr::filter(fdrset,grepl("98.0724",modified_peptide))

# Filter out protein of interest out of all butanol-modified PSMs 
#Marked_proteins <- dplyr::filter(Marked_peptides,grepl(Protein_code,protein))
Marked_proteins <- fdrset[grep(Protein_code,fdrset$protein),]

# Marked_proteins contains all butanol-modified PSMs. 
# The script now changes the score on those PSMs and condenses the multiple modified PSMs
# into a singular, modified peptide.

# Add column for protein entry name 
Sample <- data.frame(Marked_proteins$protein)

# Add column of raw E-value data.
Sample$V2 <- as.array(Marked_proteins$"e-value")

# Add column of modified PSMs. 
Sample$V3 <- as.array(Marked_proteins$modified_peptide)

# Name columns. "RPep" is the peptide sequence alone, 
names(Sample) <- c("Prot","Eval","MPep")

# Convert E-value to confidence
Sample$confidence <- as.numeric(1/(Sample$Eval))

# Remove non-butanol modifications (deamidation of N/Q and oxidation of M)
Sample$MPep <- str_replace_all(Sample$MPep,"\\[15.9949\\]","") # Remove Met ox 
Sample$MPep <- str_replace_all(Sample$MPep,"\\[0.9840\\]","") # Remove QN deamidation
Sample$MPep <- str_replace_all(Sample$MPep,"\\[57.0215\\]","") # Remove Cys alkylation, if listed as a modification
Sample$MPep <- str_replace_all(Sample$MPep,"\\[27.9950\\]","") # Remove N-formylation

# Filter out E-value more than the cutoff (less than the confidence cutoff).
Sample_CF <- subset(Sample, confidence_cutoff < confidence)

# Sum all PSMs with the same MPep string.
# Identical peptides are merged, and their confidence scores added.
Sample_Filt <- Sample_CF %>% group_by(MPep) %>% summarise_at(vars(confidence),sum)

# Clip tryptic residues off of MPep string.
# This is necessary for finding the residue number within the protein sequence.
Sample_Filt$Clip <- str_sub(Sample_Filt$MPep,3,-3)

# Create column with raw peptide (IE, no butanol labels) column. This is needed
# to match the location of the peptide to within the protein sequence and determine
# the butanol-modified residue location. 
Sample_Filt$Raw <- str_replace_all(Sample_Filt$Clip,paste("\\[",tagmass,"\\]",sep=""),"")

# Create data frame consisting of two columns.
# The "start" is the location of the butanol in the TProt sequence. The "end" sequence
# is the "start" plus the length of the detected tryptic peptide (this is useless, to my telling).
# str_locate() finds the number of the start of the "raw peptide" sequence. This is then
# added to the register of the butanol-modified residue within the peptide sequence with regexpr().
# An offset of -2 is required to keep the register-I'm not sure why.
Sample_aa_pos <- as.data.frame(str_locate(protein_sequence,Sample_Filt$Raw) + regexpr(paste("\\[",tagmass,"\\]",sep=""), Sample_Filt$Clip) + -2)

# BPeps will be the final .csv output. First, add the raw (no butanol) peptide data.
# Name the column "unmodified_sequence".
BPeps <- data.frame(Sample_Filt$Raw)
names(BPeps) <- "unmodified_sequence"

# Add modified butanol locations onto this data frame.
BPeps$modified_sequence <- Sample_Filt$Clip

# Take the butanol locations from the start location and add it onto the final data frame.
BPeps$butanol_location <- Sample_aa_pos$start

# Determine the identity of the butanolized residue, add it to a separate column (butanolized_residue)
BPeps$butanolized_residue <- substring(protein_sequence, BPeps$butanol_location,BPeps$butanol_location)

# Add "Score" value of residue.
# Technically, this isn't the confidence of each PSM-it's the summed confidence of every PSM containing
# a butanol on the residue indicated in the BPeps data frame.
# Listed as a generic "score" variable.
BPeps$score <- Sample_Filt$confidence

# Collapse identical sites.
# Some identical butanol sites occur in overlapping peptides or peptides which are not
# completely digested. This takes identical butanol locations on different peptides
# and adds them together, summing the scores.
# For example, DMSI[72.0575]LLK and DMSI[72.0575]LLKTR
# are treated as different peptides despite having the same crosslinked
# residue-this script adds the scores of them together, leaving
# one "modified" site. 
Residues <- BPeps %>% group_by(butanol_location) %>% summarize_at(vars(score),sum)

# Determine the identity of the butanolized residue, add it to a separate column (butanolized_residue)
Residues <- cbind(butanolized_residue = substring(protein_sequence, Residues$butanol_location,Residues$butanol_location), Residues)

# Write final butanol-modified residues to a csv. Can use a tab separator with "\t" instead. 
write.table(Residues, paste(OutputFileBOH,'-FDR_1of',fdrpoint,sep=""), col.names=TRUE, sep=",", row.names=FALSE)

# Write peptide data to a separate file, if desired. 
write.table(BPeps, OutputFilePeps, col.names=TRUE, sep=",", row.names=FALSE)

# Output modified residues as a single string for pymol
print(paste(Residues$butanol_location,collapse="+"))


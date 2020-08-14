#Update the fasta files especially for muting the disqualified gaps. 

library(tidyverse)
library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(dplyr)
library(stringi)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#Initial word size (Window) as 27(12+3+12). 
#convert ---AAA/AAA--- to ===AAA/AAA===
#convert ---/---       to ===/===
#convert ---AAA---     to ===AAA===


#Remove any large gaps
Filter_update= function(x,X){
  st  = start(x)
  en  = end(x)
  wid = width(x) 
  for(i in 1:length(x)){
    subseq(X,start=st[i], end=en[i]) = stri_rand_strings(1,wid[i],'[=]')
  }
  return(X)
}

#Clean up of start-stop position
start_stop_test = function(X){
  Start = substr(X,1,Window)
  IsGap.1 = grepl('-',Start)
  
  Stop = substr(X,(nchar(X)-(Window-1)),nchar(X))
  IsGap.2 = grepl('-',Stop) 
  
  if(IsGap.1 == TRUE){
    Start = gsub('-','=',Start)
    substr(X,1,Window) = Start
  }
  
  if(IsGap.2 == TRUE){
    Stop = gsub('-','=',Stop)
    substr(X,(nchar(X)-(Window-1)),nchar(X)) = Stop
  }
  return(X)
}

#Adjust gap distance
ad_gap_test = function(y,seq,ref){
  
  if(length(y)==0){return(c(seq,ref))}
  
  pos.all = start(y)
  wid.all = width(y) 
  pos.en  = end(y)
  
  for(i in 1:length(pos.all)){
    pos = pos.all[i]
    wid = wid.all[i]
    en = pos.en[i]
    left = substr(seq, start = pos-Wall-3, stop = pos-1)
    right= substr(seq, start = pos+wid, stop = pos+wid+(Wall+2))
    window= substr(ref, start = pos-Wall-3, stop = pos+wid+(Wall+2))
    
    l1 = grep("[-=]", left)
    r1 = grep("[-=]", right)
    Wid = grep("[-=]", window)
    
    if(length(l1)!=0 | length(r1)!=0 | length(Wid)!=0){
      subseq(seq,start=pos, end=en) = stri_rand_strings(1,wid,'[=]')
    }
  }
  
  return(c(seq,ref))
}


# seq1 = "AAT---AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG---ATCACCGCG---TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCATAATAGGGCCGTC---GTAATTGTCTAATATAG------ATAGTA---"
# seq2 = "TAA------AA---AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT---AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTG---CCCAATATTTAATATGATGA---TAGCTATAA"


#inFile = "../Raw_data/mapped_cds_mafft/ENSG00000000460.fa"


main = function(inFile,ouDir,num){
  #oudir="../Data/updated_cds/"
  dna = readDNAStringSet(inFile,format = "fasta")
  name = names(dna)
  
  #Variables setting
  num = as.numeric(num)
  Window  <<- num
  Wall <<- 2*Window
  
  #String mode
  spec.1 = toString(dna[[1]],width=NULL)
  spec.2 = toString(dna[[2]],width=NULL)
  
  #Find gap range
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[1]]
  r = g[[2]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  #Test gap range
  gap = c(3,6,9,12)
  m.null = m[!(wid.m %in% gap)]
  r.null = r[!(wid.r %in% gap)]
  
  m.null.l = length(m.null)
  r.null.l = length(r.null)
  
  if(m.null.l>0 & r.null.l>0){
    M = Filter_update(m.null,spec.1)
    R = Filter_update(r.null,spec.2 )  
  }else if(m.null.l>0 & r.null.l<=0){
    M = Filter_update(m.null,spec.1)  
    R = spec.2 
  }else if(m.null.l<=0 & r.null.l>0){
    M = spec.1
    R = Filter_update(r.null,spec.2) 
  }else{
    M = spec.1 
    R = spec.2 
  }
  
  
  #Clean up of start-stop area
  M.1 = start_stop_test(M)
  R.1 = start_stop_test(R)
  
  
  #Update old gap length 
  M.1.1 = str_split(as.character(M.1), '')
  R.1.1 = str_split(as.character(R.1), '')
  
  g.m = lapply(M.1.1, function(x) { IRanges(x == '-')})
  g.r = lapply(R.1.1, function(x) { IRanges(x == '-')})
  
  g.m = IRangesList(g.m)
  g.m = g.m[[1]]
  g.r = IRangesList(g.r)
  g.r = g.r[[1]]
  
  # Update the gap via swapping reference
  m.all = g.m[width(g.m)==3 | width(g.m)==6 | width(g.m)==9 | width(g.m)==12 ]
  r.all = g.r[width(g.r)==3 | width(g.r)==6 | width(g.r)==9 | width(g.r)==12 ]
  
  M_R = ad_gap_test(m.all,M.1,R.1)
  R_M = ad_gap_test(r.all,M_R[2],M_R[1])
  
  #Update the sequence
  updated_m = R_M[2]
  updated_r = R_M[1]
  
  new_seq = list(updated_m,updated_r)
  
  write.fasta(sequences = new_seq,names = name, nbchar=80,
              open = "w",as.string = TRUE, file.out = paste0(ouDir,basename(inFile)) )
  
  # write.fasta(sequences = new_seq,names = name, nbchar=80,
  #             open = "w",as.string = TRUE, file.out = "Test_update/test_output.fa" )
  
}

#compatible with run_test.R
#if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])
#}

###############################################


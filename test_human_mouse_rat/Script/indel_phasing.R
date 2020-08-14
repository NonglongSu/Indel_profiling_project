#Calculate the freq distribution of phase 0 / phase 1 / phase 2 INDELs. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")
#file="Test_pseudo_seq/tmp/test2.fa"

phasing = function(file){
  
  dna = readBStringSet(file,format = "fasta")
  
  #String mode
  M = toString(dna[[1]],width=NULL)
  R = toString(dna[[2]],width=NULL)
  
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[1]]
  r = g[[2]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  
  df.m = data.frame("pos"=pos.m,"wid"=wid.m)
  df.r = data.frame("pos"=pos.r,"wid"=wid.r)
  
  if(l.m > 0 & l.r == 0){
    return(df.m)
  }else if(l.m == 0 & l.r > 0){
    return(df.r)
  }else if(l.m > 0 & l.r > 0){
    df = merge(df.m,df.r,all=TRUE)
    return(df)
  }else{
    return(NULL)
  }
  
}


Record = function(PST){
  if(is.null(PST)){return(PST)}
  rem   = PST%%3
  phase.0 = length(which(rem==1))
  phase.1 = length(which(rem==2)) 
  phase.2 = length(which(rem==0)) 
  
  df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
  return(df.phase)
}  

# file1 = "../Data/mapped_cds_Out/ENSG00000049246.fa"
# file1 = "../Data/mapped_cds_Out/ENSG00000000460.fa"


#m1 = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG------ATCACCGCGAGATCATGTCTAGTTATGAACGGC======GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG======ATAGTA"
#r1 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"
# 
# m2 = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCATAATAG---GGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
# r2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGATC---GCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

main=function(dir,ouFile,ouFig){
  
  File.total = list.files(dir, full.names = TRUE)
  
  PST = c()
  for(i in 1:length(File.total)){
    phase_score = phasing(File.total[i])
    PST =rbind(PST,phase_score)
  }
  
  
  pos.3 = PST[PST$wid==3,]$pos
  pos.6 = PST[PST$wid==6,]$pos
  pos.9 = PST[PST$wid==9,]$pos
  pos.12 = PST[PST$wid==12,]$pos
  
  phase.3  = Record(pos.3)
  phase.6  = Record(pos.6)
  phase.9  = Record(pos.9)
  phase.12 = Record(pos.12)
  
  Phase.df = rbind(phase.3,phase.6,phase.9,phase.12)
  Phase.df = rbind(Phase.df,colSums(Phase.df))
  gap_length = c('3','6','9','12',"T")
  Phase.df = cbind(Phase.df,gap_length)
  
  write.table(Phase.df,file = ouFile,sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  #dev.off(dev.list())
  #Pie plot
  #outFig="../Figure/Phase_Pie_chart.pdf"
  pdf(ouFig)
  
  name    = colnames(Phase.df)[1:3]
  p.3 = as.numeric(Phase.df[1,1:3])
  p.6 = as.numeric(Phase.df[2,1:3])
  p.9 = as.numeric(Phase.df[3,1:3])
  p.12 = as.numeric(Phase.df[4,1:3])
  p.T = as.numeric(Phase.df[5,1:3])
  
  par(mfrow=c(3,2))
  
  #Sometimes if all values in P._ is zero, then errors shows.  
  pie(p.3,labels=name,main = "Pie chart of Three phases of gap length 3",col=rainbow(length(p.3)))
  pie(p.6,labels=name,main = "Pie chart of Three phases of gap length 6",col=rainbow(length(p.6)))
  pie(p.9,labels=name,main = "Pie chart of Three phases of gap length 9",col=rainbow(length(p.9)))
  pie(p.12,labels=name,main = "Pie chart of Three phases of gap length 12",col=rainbow(length(p.12)))
  pie(p.T,labels=name,main = "Pie chart of Three phases of gap length of all",col=rainbow(length(p.T)))
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])

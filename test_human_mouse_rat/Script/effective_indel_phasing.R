#Calculate the freq distribution of effective phase 1 & phase 2 INDELs. 
#Define effective as "can cause a substitution".

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)
library(dplyr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")
#file="Test_pseudo_seq/tmp/test2.fa"

eff_phase = function(seq1,seq2,pos,wid,len){#seq1--reference
  seq2.char = str_split(seq2,"")[[1]]
  eff.sub.1 = 0
  non.sub.1 = 0
  eff.sub.2 = 0
  non.sub.2 = 0
  for(j in 1: len){
    if(pos[j] %%3 == 1){#phase-0
      next
    }else if(pos[j] %%3 == 0){#phase-2
      pos.ori = pos[j] - 2 
      unit.1 = substr(seq1,pos.ori,pos.ori+2)
      unit.2 = substr(seq1,pos.ori+wid[j],pos.ori+wid[j]+2)
      sub = paste0(c(seq2.char[pos[j]-2],seq2.char[pos[j]-1],seq2.char[pos[j]+wid[j]]),collapse = "")
      sec = codon[[which(sapply(codon,function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        non.sub.2 = non.sub.2 + 1
      }else{
        eff.sub.2 = eff.sub.2 + 1
      }
    }else{#phase-1
      pos.ori = pos - 1 
      unit.1 = substr(seq1,pos.ori,pos.ori+2)
      unit.2 = substr(seq1,pos.ori+wid[j],pos.ori+wid[j]+2)
      sub = paste0(c(seq2.char[pos[j]-1],seq2.char[pos[j]+wid[j]],seq2.char[pos[j]+wid[j]+1]),collapse = "")
      sec = codon[[which(sapply(codon,function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        non.sub.1 = non.sub.1 + 1
      }else{
        eff.sub.1 = eff.sub.1 + 1
      }
    }
  }  
 res = c(eff.sub.1,non.sub.1,eff.sub.2,non.sub.2) 
 return(res)
}

# file1 = "../Data/mapped_cds_Out/ENSG00000049246.fa"
# file1 = "../Data/mapped_cds_Out/ENSG00000000460.fa"

# dir = "~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data/mafft/mapped_cds_Out"

# m1 = "AAT---AAAC---AAAGAAT---G"
# r1 = "TAAATCATTAAAAAGAATTTGATG"
# 
# m2 = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCATAATAG---GGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
# r2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGATC---GCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

main=function(dir,ouFile,ouFig){
  #Create a codon table
  codon <<- list (c("TTT","TTC"),
                  c("TTA","TTG","CTT","CTC","CTA","CTG"),
                  c("ATT","ATC","ATA"),
                  c("ATG"),
                  c("GTT","GTC","GTA","GTG"),
                  c("TCT","TCC","TCA","TCG","AGT","AGC"),
                  c("CCT","CCC","CCA","CCG"),
                  c("ACT","ACC","ACA","ACG"),
                  c("GCT","GCC","GCA","GCG"),
                  c("TAT","TAC"),
                  c("CAT","CAC"),
                  c("CAA","CAG"),
                  c("AAT","AAC"),
                  c("AAA","AAG"),
                  c("GAT","GAC"),
                  c("GAA","GAG"),
                  c("TGT","TGC"),
                  c("TGG"),
                  c("CGT","CGC","CGA","CGG","AGA","AGG"),
                  c("GGT","GGC","GGA","GGG"),
                  c("TAA","TGA","TAG")
  )
  
  M.eff.1 = 0 ; R.eff.1 = 0
  M.non.1 = 0 ; R.non.1 = 0
  M.eff.2 = 0 ; R.eff.2 = 0
  M.non.2 = 0 ; R.non.2 = 0
  
  File.total = list.files(dir, full.names = TRUE)
  for(i in 1:length(File.total)){
    dna = readBStringSet(File.total[i],format = "fasta")
    
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
    
    if(l.m > 0){
      eff_M = eff_phase(R,M,pos.m,wid.m,l.m)
      M.eff.1 = M.eff.1 + eff_M[1]
      M.non.1 = M.non.1 + eff_M[2]
      M.eff.2 = M.eff.2 + eff_M[3] 
      M.non.2 = M.non.2 + eff_M[4] 
    }
    if(l.r > 0){
      eff_R = eff_phase(M,R,pos.r,wid.r,l.r)
      R.eff.1 = R.eff.1 + eff_R[1]
      R.non.1 = R.non.1 + eff_R[2]
      R.eff.2 = R.eff.2 + eff_R[3] 
      R.non.2 = R.non.2 + eff_R[4] 
    }
  }
  
  #Generate an effective phase table
  focal_name = c('Mouse','Rat')
  Phase.M = data.frame("eff_1"=M.eff.1,"non_1"=M.non.1,"eff_2"=M.eff.2,"non_2"=M.non.2)
  Phase.R = data.frame("eff_1"=R.eff.1,"non_1"=R.non.1,"eff_2"=R.eff.2,"non_2"=R.non.2)
  Phase.df= rbind(Phase.M,Phase.R)
  Pdf = cbind(focal_name, Phase.df)
  
  write.table(Pdf,file = ouFile,sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  
  #dev.off(dev.list())
  #Pie plot
  #ouFig="../Figure/Phase_Pie_chart.pdf"
  
  P.M.1 = t(select(Phase.M,eff_1,non_1))
  P.M.2 = t(select(Phase.M,eff_2,non_2))
  P.R.1 = t(select(Phase.R,eff_1,non_1))
  P.R.2 = t(select(Phase.R,eff_2,non_2))
  col.name = c("phase_1","phase_2")
  
  pdf(file=ouFig)
  par(mfrow=c(2,1))
  barplot(cbind(P.M.1,P.M.2),main = "effective phase distribution of mouse", names.arg = col.name, col = rainbow(2),bty='L' )
  legend(2.2,0,legend = c("effective","non-effective"),xpd=TRUE, cex=0.4, fill = rainbow(2))
  barplot(cbind(P.R.1,P.R.2),main = "effective phase distribution of rat", names.arg = col.name, col = rainbow(2))
  legend(2.2,0,legend = c("effective","non-effective"),xpd=TRUE, cex=0.4, fill = rainbow(2))
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])

#Generate the Edge bias source file(txt)

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(stringr)
library(tidyr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

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
  
  pos.m = start(m)
  pos.r = start(r)
  
  df.m = data.frame("pos"=pos.m,"wid"=wid.m,"file"=rep(file,length(pos.m)))
  df.r = data.frame("pos"=pos.r,"wid"=wid.r,"file"=rep(file,length(pos.r)))
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
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

Extract = function(file,idx,wid){
  dna   = readBStringSet(file,format = "fasta")
  x = toString(dna[[1]],width=NULL)
  y = toString(dna[[2]],width=NULL)
  
  start = idx-Window
  stop  = idx+Window+wid-1
  wid.1 = substr(x, start, stop)
  wid.2 = substr(y, start, stop)
  
  gap = grep('-',wid.1)
  gap.l = length(gap)
  
  if(gap.l>0){
    swap(wid.1,wid.2)
  }
  df = data.frame("wid.1"=wid.1,"wid.2"=wid.2)
  return(df)
}


table_magic = function(x,k,y,z){
  pos.edge  = x$pos[which(z==y)]
  wid.edge  = x$wid[which(z==y)]
  
  file.list = k$file[which(z==y)]
  file.list = as.character(file.list)
  
  windows = c()
  for(i in 1:length(file.list)){
    val  = Extract(file.list[i],pos.edge[i],wid.edge[i]) 
    windows = rbind(windows,val)
  }
  
  amino_acid.1 = Biostrings::translate(DNAStringSet(windows$wid.1),genetic.code=GENETIC_CODE,if.fuzzy.codon="solve")
  windows$amino_acid.1 = as.character(amino_acid.1)
  
  tmp.1 = str_split(windows$wid.2,"")
  amino_acid.2 = list()
  for (i in 1:length(tmp.1)){
    tmp.2 = translate(tmp.1[[i]],NAstring = "-", ambiguous = FALSE)
    tmp.2 = paste(tmp.2, collapse=" ")
    amino_acid.2[i] = tmp.2
  }
  
  amino_acid.2 = as.character(amino_acid.2)
  windows$amino_acid.2 = amino_acid.2
  
  return(windows)
}

#dir1 = "../Data/mapped_cds_plus"
#dir2 = "../Data/mapped_cds_Out"

main=function(dir1,dir2,ouDir,num){
  
  File1.total = list.files(dir1, full.names = TRUE)
  File2.total = list.files(dir2, full.names = TRUE)
  
  #mapped_plus
  PST.1 = c()
  for(i in 1:length(File1.total)){
    pScore.1 = phasing(File1.total[i])
    if(!is.null(pScore.1)){
      PST.1 = rbind(PST.1,pScore.1)
    }
  }
  #mapped_cds_Out
  PST.2 = c()
  for(i in 1:length(File2.total)){
    pScore.2 = phasing(File2.total[i])
    if(!is.null(pScore.2)){
      PST.2 =rbind(PST.2,pScore.2)
    }
  }
  
  #Perfecto: The length of PST.1 equal to length of PST.2. 
  
  #set up edge
  num = as.numeric(num) 
  Window <<- num 
  Edge = c(-Window,Window)                           
  diff.T = PST.2$pos - PST.1$pos 
  
  
  #left bias
  align_ori.l    = table_magic(PST.1,PST.1,Edge[1],diff.T) 
  align_better.l = table_magic(PST.1,PST.2,Edge[1],diff.T)  
  #right bias
  align_ori.r    = table_magic(PST.1,PST.1,Edge[2],diff.T) 
  align_better.r = table_magic(PST.1,PST.2,Edge[2],diff.T)  
  
  # align.l = data.frame("AA_ref"=align_ori.l$amino_acid.1,"AA_ori_align"=align_ori.l$amino_acid.2,"AA_better_align"=align_better.l$amino_acid.2,
  #                      "DNA_ref"=align_ori.l$wid.1,"DNA_ori_align" = align_ori.l$wid.2,"DNA_better_align"=align_better.l$wid.2)
  # 
  # align.r = data.frame("AA_ref"=align_ori.r$amino_acid.1,"AA_ori_align"=align_ori.r$amino_acid.2,"AA_better_align"=align_better.r$amino_acid.2,
  #                      "DNA_ref"=align_ori.r$wid.1,"DNA_ori_align" = align_ori.r$wid.2,"DNA_better_align"=align_better.r$wid.2)
  
  #Write out
  ouFile.1 = "align.ori.l.txt"
  ouFile.2 = "align.better.l.txt"
  ouFile.3 = "align.ori.r.txt"
  ouFile.4 = "align.better.r.txt"
 
  # ouFile.1 = "align.l.txt" 
  # ouFile.2 = "align.r.txt" 
  
  align.edge = list(align_ori.l,align_better.l,align_ori.r,align_better.r)
  ouFile     = c(ouFile.1,ouFile.2,ouFile.3,ouFile.4)
  
  for(i in 1:length(ouFile)){
    write.table(align.edge[i], file = paste0(ouDir,basename(ouFile[i])),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  }
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3],args[4])



#wind.l = unite(windows.l,"window",c(amino_acid.1,wid.1,wid.2,amino_acid.2),sep='\n',remove = TRUE)
#windows.l$wid.12 = paste(windows.l$wid.1,windows.l$wid.2, sep = '\n')

#Generate a list of files that "not_multiple of 3" and "species inside less than 3". 

library(Biostrings)
library(scales)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#dir = "../Raw_data.2.outgroups/cds_seq"
#output = "../Raw_data/summary_perc.txt"

cal_perc = function(dir,num,ouFile.1,ouFile.2){
  
  num = as.numeric(num)
  L  <<- num
  
  allData = list.files(dir,full=TRUE)
  for(i in allData){
    cds = readDNAStringSet(i,format="fasta")
    width = width(cds)
    perc = sum(width %%3==0)/length(cds)
    
    if(length(cds) < L){
        cat(as.character(basename(i)),file = ouFile.1,sep = "\n",append = TRUE)
    }else{
      if(perc < 1){
        cat(as.character(basename(i)),file = ouFile.2,sep = "\n",append = TRUE)
      }
    }
  }
}

args = commandArgs(trailingOnly = TRUE)
cal_perc(args[1],args[2],args[3],args[4])

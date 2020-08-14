library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

inDir  = "../Tmp1"

test_early_stop = function(inDir,ouFile,spec){
  
  File.holder = list.files(inDir,full=TRUE)
  #spec = c("human","hamster","mouse","rat")
  
  
  stop.codon = c("TAG","TGA","TAA")
  fileList   = c()
  filename   = c()
  for( i in File.holder){
    dna = readDNAStringSet(i,format="fasta")
    len = width(dna)
    for(j in 1:length(dna)){
      pos   = gregexpr(paste(stop.codon,collapse="|"),dna[[j]])[[1]]
      pos.1 = pos[pos %% 3 == 1]
        if(length(pos.1) > 0 && pos.1[1]<(len[j]-2) ){#Always targetting the 1st stop codon
              #dna[[j]]  = dna[[j]][-(pos.1[1]+3):-(width(dna)[j])]
              #The stop codon can be rare codon
              filename  = basename(i)
              fileList = c(fileList,filename)
              break
            }
    }
  }
  
    #ouFile = "../Raw_data.2.outgroups/early_stopCodon.txt"
    cat(fileList,file = ouFile,sep = "\n",append = FALSE)
}


args = commandArgs(trailingOnly = TRUE)
test_early_stop(args[1],args[2],eval(parse(text=args[3])))



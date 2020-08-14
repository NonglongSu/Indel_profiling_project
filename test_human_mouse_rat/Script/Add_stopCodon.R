#Add stopcodon to the end of the sequence who does not have it. 

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")
library(Biostrings)
library(stringr)
library(seqinr)

#file_name = "../Raw_data.2.outgroups/cds_seq/ENSG00000084072.fa"
#output = "../Raw_data.2.outgroups/cds_seq_stop/"

stopCodon_add = function(file_name,output){
  cds = readDNAStringSet(file_name,format="fasta")
  wid = width(cds)
  
  new_seq = list()
  
  for(i in 1:length(cds)){
    old = toString(cds[[i]],width=NULL)
    ending = wid[i]
    Stop = str_sub(old,start=wid[i]-2,ending)
    if(Stop != "TAG" & Stop != "TGA" & Stop != "TAA"){
      old = gsub('^(.*)$', '\\1---\\2', old)
    }
    new_seq[i]=old
  }

name = names(cds)
write.fasta(sequences = new_seq,names = name, nbchar=80,
              open = "w",as.string = TRUE, file.out = paste0(output,basename(file_name)))

}

args = commandArgs(trailingOnly = TRUE)
stopCodon_add(args[1],args[2])


############################
#file_name = "../Raw_data/haha.fa"
# cds <- readDNAStringSet(inFile,format="fasta")


#Abnormal alignment: align amino acid to DNA.

library(seqinr)

## The file within the folder order lowest to sort. 
#file1 = "mafft_out/ENSG00000000460.fa"
#file2 = "cds_seq_stop/ENSG00000000460.fa"

map_cds2aa = function(file1,file2,file3){

#file1 = file.path("..","Raw_data",file1)
#file2 = file.path("..","Raw_data",file2)

#file1 = file.path("..","Raw_data.2.outgroups",file1)
#file2 = file.path("..","Raw_data.2.outgroups",file2)


#print(c(file1,file2))
data1 = read.fasta(file1,seqtype = "AA",as.string = FALSE,
                   set.attributes = FALSE)
data2 = read.fasta(file2,seqtype = "DNA",as.string = FALSE,
                   set.attributes = FALSE,forceDNAtolower = FALSE)

for(i in 1:length(data1)){
  if(length(which(data1[[i]]=="-")) != 0){
    pos = which(data1[[i]]=="-")
    for(j in 1: length(pos)){
      data2[[i]] = append(data2[[i]],c("-","-","-"),after=(pos[j]-1)*3)
    }
  }
}



name = names(data1)


#file3 = "../Raw_data/mapped_cds_seq/"
write.fasta(sequences = data2,names = name, nbchar=80,
            open = "w",as.string = FALSE, file.out = paste0(file3,basename(file1)))
}

args = commandArgs(trailingOnly = TRUE)
map_cds2aa(args[1],args[2],args[3])

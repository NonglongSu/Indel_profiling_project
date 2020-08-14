library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script/QC")

inDir = "../Raw_data.2.outgroups/cds_seq"

test_N = function(inDir,ouFile){
    File.holder = list.files(inDir,full=TRUE)
    
    N.list = c()

    for(i in 1: length(File.holder)){
        dna   = readDNAStringSet(File.holder[i],format="fasta")
        len   = length(grep('N', dna, ignore.case = TRUE))
        if(len > 0){
           N.list = c(N.list,basename(File.holder[i]))
        }
    }

    cat(N.list,file = ouFile,sep = "\n",append = FALSE)
}

args = commandArgs(trailingOnly=TRUE)
test_N(args[1],args[2])

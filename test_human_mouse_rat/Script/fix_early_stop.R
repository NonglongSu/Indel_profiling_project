#Replace the stop codon symbol(*) with X. 

library(Biostrings)
library(stringr)
library(seqinr)
#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

# file_name = "../Raw_data.2.outgroups/aa_seq/ENSG00000073169.fa"
# file.out  = "../Raw_data.2.outgroups/aa_seq_fix/"

stopCodon_fix = function(file_name,file.out){
  aa = readAAStringSet(file_name,format="fasta")
  
  tree.new = list()
  for (i in 1:length(aa)) {
    tree.new[[i]] = str_replace_all(toString(aa[[i]]), "\\*", "X")
  }

  write.fasta(tree.new, names(aa), paste0(file.out,basename(file_name)), open = "w", nbchar = 80, as.string = FALSE)
}

args = commandArgs(trailingOnly = TRUE)
stopCodon_fix(args[1],args[2])

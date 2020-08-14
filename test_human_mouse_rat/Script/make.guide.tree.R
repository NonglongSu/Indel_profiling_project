library(ape)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#Url keeps updating
url = "ftp://ftp.ensembl.org/pub/current_compara/species_trees/38_mammals_EPO_default.nh"

#keep = c("homo_sapiens", "mus_musculus_strain_reference_cl57bl6", "rattus_norvegicus", "cricetulus_griseus")
 keep = c("homo_sapiens","mus_musculus_strain_reference_cl57bl6", "rattus_norvegicus")

a = read.tree(url)

a$tip.label = str_to_lower(a$tip.label) 

drop = !(a$tip.label %in% keep)

b = drop.tip(a, a$tip.label[drop])

#b$tip.label = c("hamster","mouse","rat","human")
 b$tip.label = c("mouse","rat","human")


cat(write.tree(b))
 
 
 
 

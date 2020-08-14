#extract the gene ID from ensembl database

# rm(list=ls())
# invisible(dev.off())
######################
# .libPaths()
# [1] "/home/zzhu49/R/x86_64-pc-linux-gnu-library/3.4"
# [2] "/usr/local/lib/R/site-library"                 
# [3] "/usr/lib/R/site-library"                       
# [4] "/usr/lib/R/library"  

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocInstaller")
# biocLite("Biostrings")
# biocLite("biomaRt")
# biocLite("S4Vectors")
# biocLite("BiocGenerics")
# install.packages("readr")
# install.packages("dplyr")

library(BiocInstaller)
library(BiocGenerics)
library(S4Vectors)
library(Biostrings)
library(biomaRt)
library(readr)
library(dplyr)

output = "../Raw_data.2.outgroups/geneId.txt"

ensembl=useMart("ensembl")
human = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# b = listAttributes(human)
# b$name
# a = listFilters(human)
# a$name

attributes = c("ensembl_gene_id",
               "cchok1gshd_homolog_ensembl_gene",
               "cchok1gshd_homolog_orthology_type",
               "mmusculus_homolog_ensembl_gene",
               "mmusculus_homolog_orthology_type" ,
               "rnorvegicus_homolog_ensembl_gene",
               "rnorvegicus_homolog_orthology_type"
)


filters= "with_ccds"


orth = getBM(attributes,filters,values=TRUE,
             mart = human, uniqueRows=TRUE)

#  element-wise comparison single &. 
orth_filtered = orth %>% filter(`mmusculus_homolog_orthology_type`=="ortholog_one2one" & 
                                  `rnorvegicus_homolog_orthology_type`=="ortholog_one2one" &
                                  `cchok1gshd_homolog_orthology_type`=="ortholog_one2one")

orth_com = orth_filtered[,c(1,2,4,6)]

write.table(orth_com,output,sep='\t',row.names=F)




#Extract the CDS from ensembl API.

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)


server = "http://rest.ensembl.org"

output   = "../Raw_data/cds_seq/"                         #make cds
output.1 = "../Raw_data.2.outgroups/cds_seq/Tmp1/"        #make cds.1
output.2 = "../Raw_data.2.outgroups/cds_seq/Tmp2/"        #make cds.2

#name="ENSG00000158796"
main = function(name){
      name = gsub('\"', "", name, fixed = TRUE)
      ext  = paste("/lookup/id/",name,"?expand=1",sep = "")
      r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
      stop_for_status(r)
      id   = fromJSON(toJSON(content(r)))
      id_  = id$Transcript[["id"]][id$Transcript$is_canonical==1]
      #print(id)
      
      ext_1 = paste("/sequence/id/",id_,"?type=cds",sep = "")
      r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
      stop_for_status(r_1)
      cds_  = (content(r_1))
      #cds_ =noquote(strsplit(cds, split = "\n", fixed = TRUE))
      #print (cds_)
      
      list_name = c(name,cds_)
      return (list_name)
}

args = commandArgs(trailingOnly = TRUE)
cds.df = data.frame(x = numeric(length(args)),row.names = NULL,check.rows = FALSE, stringsAsFactors = FALSE)

for(i in 1:length(args)){
  list.name = main(args[i])
  if(i==1){
     cds.name = list.name[1]
  }
  cds.seq = list.name[2]
  cds.df$x[i] = cds.seq
}

write.table(cds.df$x,paste0(output,cds.name,".fa"),quote = FALSE,
                row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\n")
  



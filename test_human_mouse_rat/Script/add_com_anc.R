library(Biostrings)
library(stringr)
library(seqinr)

# setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#Where the outgroups provide no information will be removed 
#and this can lead to the underestimation of nucleotide variation.

#convert string to a vector.
Str2Vec = function(s){
  s1 = str_split(s,"")
  s1 = s1[[1]]
  return(s1)
}

#Dealing with codon without gaps
Find_com_anc_good = function(seq1,seq2,Anc,site,Tag){     
  init = Tag #tag are not supposed to increase unless we get legal codon set. 
  for(w in 1:length(site)){#loop of 3.
    most = max(as.data.frame(table(site[[w]]))$Freq) #Find how many times the most frequent base occurs
    diff = length(unique(site[[w]]))                 #Find how many different bases in this tree-set 
    df   = as.data.frame(table(site[[w]])) 
    var  = as.character(df[[1]])                     #Find all unique bases
      
    if(any(site[[w]] %in% '=') || any(site[[w]] %in% '+')  || most==1 ){# A/T/G/C , =/A/A/A , +/T/T/T
        seq1 = seq1[-init:(-init-2)]    #remove the unsatisfied codon
        seq2 = seq2[-init:(-init-2)]
        Anc  = Anc [-init:(-init-2)]
        Tag  = init                      
        break
      }else{
        if(most==2 && diff==2){
          if(site[[w]][1]!=site[[w]][2]){# A/T/A/T , A/T/T/A
             seq1 = seq1[-init:(-init-2)]
             seq2 = seq2[-init:(-init-2)]
             Anc  = Anc [-init:(-init-2)]
             Tag  = init 
             break 
          }else{# A/A/G/G
             base = site[[w]][1]
             Anc  = c(Anc,base)
             Tag  = Tag + 1 
          }
        }else{# A/A/G/T , A/A/A/G, A/A/A/A
          base = var[which(df$Freq == max(df$Freq))]  #pick the most frequent base as the ancestor base
          Anc  = c(Anc,base) 
          Tag  = Tag + 1 
        }
      }
    }
      
  res = list(seq1,seq2,Anc,Tag)
  return(res)
}

#Dealing with codon containing gaps
Find_com_anc_bad = function(seq1,seq2,Anc,site,Tag,Wind){
  init = Tag
  for(w in 1:Wind){
    most = max(as.data.frame(table(site[[w]]))$Freq)
    diff = length(unique(site[[w]]))
    df   = as.data.frame(table(site[[w]]))
    var  = as.character(df[[1]])
    if(most==2 && diff==2){
      if(site[[w]][1]!=site[[w]][2]){# A/T/A/T , A/T/T/A , A/-/A/- , A/-/-/A
        seq1 = seq1[-init:-(init+Wind-1)]     
        seq2 = seq2[-init:-(init+Wind-1)]     
        Anc  = Anc [-init:-(init+Wind-1)]     
        Tag  = init
        break
      }else{# A/A/G/G , A/A/-/- , -/-/A/A 
        base = site[[w]][1]
        Anc  = c(Anc,base)
        Tag  = Tag + 1
      }
    }else if(most == 1 ){# A/G/C/T, A/-/G/T
      seq1 = seq1[-init:-(init+Wind-1)]     
      seq2 = seq2[-init:-(init+Wind-1)]   
      Anc  = Anc [-init:-(init+Wind-1)]   
      Tag  = init
      break
    }else{# A/A/A/-, -/-/A/-, A/T/-/A, -/-/A/T, A/A/A/A , 
        base = var[which(df$Freq == max(df$Freq))]
        Anc  = c(Anc,base) 
        Tag  = Tag + 1
      }
    }
      
  res = list(seq1,seq2,Anc,Tag)
  return(res)
}
 

# inFile = "../Data_6.2/Mafft_Mul/mapped_cds_out/ENSG00000111319.fa"
# inFile = "../Data_6.2/sample.fa"
# ouDir  = "../Tmp2/"
# H = "AA---A"
# C = "AA---A"
# M = "TT---T"
# R = "ATGATA"

main = function(inFile,ouDir){
  
  dna = readBStringSet(inFile,format = "fasta")
  #Use the minimum length because the length might not be equal(fix in the future).
  len = min(width(dna))
  
  H = toString(dna[[1]],width=NULL)
  C = toString(dna[[2]],width=NULL)
  M = toString(dna[[3]],width=NULL)
  R = toString(dna[[4]],width=NULL)
  
  u = Str2Vec(H)
  v = Str2Vec(C)
  p = Str2Vec(M)
  q = Str2Vec(R)
  Tree=list(u,v,p,q)
  
  #Grap all legit gaps
  gaps      = lapply(Tree,function(x){IRanges(x == '-')})
  freeGaps  = unlist(IRangesList(gaps))
  
  #Find commomn ancestor and update the focal sequences. 
  p1  = p
  q1  = q
  tag = i = 1 
  anc = c()
  while(i < len){
    cSite = list()
    #Start with the first codon set
    for(j in i:(i+2)){
      Site = unlist(lapply(Tree, `[[`,j),use.names=TRUE)
      cSite = c(cSite,list(Site))
    }
    #Create a decision-maker
    Bell = any(lapply(cSite, function(x){any(x  %in% '-')}) == TRUE)
    if(Bell){
      for (k in 1:3) {#Loop of 3
        if(any(cSite[[k]] == '-')){#Determine the correct gap position
          gap.pos = i+k-1
          gap.wid = width(freeGaps[which(start(freeGaps) == gap.pos)])[1]
          if(gap.pos %% 3 == 1){#phase-0
             Window = i:(i+gap.wid-1) #between-codon gaps
          }else{
             Window = i:(i+gap.wid+2) #inside-codon gaps + edges
          }
          break
        }
      }
      #Apply the window 
      dSite = list()
      for(m in Window){
            Site  = unlist(lapply(Tree, `[[`,m),use.names=TRUE)
            dSite = c(dSite,list(Site))
      }
      res = Find_com_anc_bad(p1,q1,anc,dSite,tag,length(Window))
      i   = tail(Window,1) + 1
    }else{
      res = Find_com_anc_good(p1,q1,anc,cSite,tag)
      i   = i + 3
    }
    p1  = res[[1]]
    q1  = res[[2]]
    anc = res[[3]]
    tag = res[[4]] 
}
 
  p2   = paste(p1,collapse ="")
  q2   = paste(q1,collapse ="")
  anc2 = paste(anc,collapse="")
  
  newTree = list(anc2,p2,q2)
  
  name = c("ancestor","mouse","rat")
  write.fasta(sequences=newTree,names=name,nbchar=80,
              open ="w",as.string=TRUE,file.out=paste0(ouDir,basename(inFile)))
}

args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2])


  
  
 
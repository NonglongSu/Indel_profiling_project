#Calculate the Displacement distribution of INDELs of different lengths. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

phasing = function(file){
  
  dna = readBStringSet(file,format = "fasta")
  #String mode
  M = toString(dna[[1]],width=NULL)
  R = toString(dna[[2]],width=NULL)
  
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[1]]
  r = g[[2]]
  
  wid.m = width(m)
  wid.r = width(r)

  pos.m = start(m)
  pos.r = start(r)

  df.m = data.frame("pos"=pos.m,"wid"=wid.m)
  df.r = data.frame("pos"=pos.r,"wid"=wid.r)

  l.m = length(wid.m)
  l.r = length(wid.r)

  if(l.m > 0 & l.r == 0){
    return(df.m)
  }else if(l.m == 0 & l.r > 0){
    return(df.r)
  }else if(l.m > 0 & l.r > 0){
    df = merge(df.m,df.r,all=TRUE)
    return(df)
  }else{
    return(NULL)
  }

  # wid.T = c(wid.m,wid.r)
  # flag = wid.T %in% c(3,6,9,12)
  # if(all(flag)){
  #   return(NULL)
  # }else{
  #   return(file)
  # }
}

#dir1 = "../Data/mapped_cds_plus"
#dir2 = "../Data/mapped_cds_Out"

main=function(dir1,dir2,ouFile,ouFig,num){
  
  File1.total = list.files(dir1, full.names = TRUE)
  File2.total = list.files(dir2, full.names = TRUE)
  
  # Ab_file = c()
  # for(i in 1:length(File2.total)){
  #   filename = phasing(File2.total[i])
  #   Ab_file =rbind(Ab_file,filename)
  # }
  
  # wid.1 = PST.1$wid
  # wid.1[!(wid.1 %in% c(3,6,9,12))]
  # 
  # wid.2 = PST.2$wid
  # PST.1$pos[!(wid.2 %in% c(3,6,9,12))]
  
  PST.1 = c()
  for(i in 1:length(File1.total)){
    pScore.1 = phasing(File1.total[i])
    PST.1 =rbind(PST.1,pScore.1)
  }
  
  PST.2 = c()
  for(i in 1:length(File2.total)){
    pScore.2 = phasing(File2.total[i])
    PST.2 =rbind(PST.2,pScore.2)
  }
  
  #Perfecto: The length of PST.1 equal to length of PST.2. 
  #mapped_cds_plus
  pos1.3 =  PST.1[PST.1$wid==3,]$pos
  pos1.6 =  PST.1[PST.1$wid==6,]$pos
  pos1.9 =  PST.1[PST.1$wid==9,]$pos
  pos1.12 = PST.1[PST.1$wid==12,]$pos
  #mapped_cds_Out
  pos2.3 =  PST.2[PST.2$wid==3,]$pos
  pos2.6 =  PST.2[PST.2$wid==6,]$pos
  pos2.9 =  PST.2[PST.2$wid==9,]$pos
  pos2.12 = PST.2[PST.2$wid==12,]$pos
  
  diff.3  = pos2.3  - pos1.3
  diff.6  = pos2.6  - pos1.6
  diff.9  = pos2.9  - pos1.9
  diff.12 = pos2.12 - pos1.12
  
  #Create diplacement column
  window = as.numeric(num)
  
  Dis = seq(-window,window,1)
  G.3 = c()
  G.6 = c()
  G.9 = c()
  G.12 = c()
  
  for(i in 1:length(Dis)){
    G.3[i]  = length(diff.3[diff.3==Dis[i]])
    G.6[i]  = length(diff.6[diff.6==Dis[i]])
    G.9[i]  = length(diff.9[diff.9==Dis[i]])
    G.12[i] = length(diff.12[diff.12==Dis[i]])
  }
  GD.df = rbind(G.3,G.6,G.9,G.12)
  colnames(GD.df) = Dis
  rownames(GD.df) = c(3,6,9,12)
  
  T = c()
  for(i in 1:length(Dis)){
    T[i] = sum(GD.df[,i])
  }
  GD.df= rbind(T,GD.df)
  
  # #Absolute displacement
  # Dis.abs = c(3,2,1,0)
  # GD.abs=c()
  # for(j in 1:ceiling(length(Dis)/2) ){
  #   if(j != length(Dis)-j+1){
  #     GD.abs = rbind(GD.abs,GD.df[,j] + GD.df[,length(Dis)-j+1])
  #   }else{
  #     GD.abs = rbind(GD.abs,GD.df[,j])
  #   }
  # }
  # rownames(GD.abs) = Dis.abs
  
  #Write out files
  #ouFile = "../Data/Displacement_Power.txt"
  write.table(GD.df,file = ouFile,sep = "\t",append = FALSE,quote = FALSE,row.names = TRUE,col.names = NA)
  
  
  
  #dev.off(dev.list())
  #Bar chart
  #outFig="../Figure/Displacement_Power_PieChart.pdf"
  pdf(file = ouFig)
  
  lengths    = rownames(GD.df)
  displacement = colnames(GD.df)
  GD = t(GD.df)
  barplot(GD,main = "Displacement distribution of different gaps", names.arg = lengths, xlab = "Lengths", ylab = "Displacement",
          col = rainbow(length(displacement)))
  
  legend("topright",displacement,cex=1,fill = rainbow(length(displacement)))
  
  # par(mfrow=c(3,2))
  # pie(GD.abs[,1],labels=name,main  = "Pie chart of displacement of gap length 3",col=rainbow(length(GD.abs[,1])))
  # pie(GD.abs[,2],labels=name,main  = "Pie chart of displacement of gap length 6",col=rainbow(length(GD.abs[,2])))
  # pie(GD.abs[,3],labels=name,main  = "Pie chart of displacement of gap length 9",col=rainbow(length(GD.abs[,3])))
  # pie(GD.abs[,4],labels=name,main = "Pie chart of displacement of gap length 12",col=rainbow(length(GD.abs[,4])))
  # pie(GD.abs[,5],labels=name,main  = "Pie chart of displacement of gap length of all",col=rainbow(length(GD.abs[,5])))
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3],args[4],args[5])

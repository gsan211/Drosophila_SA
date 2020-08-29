#!/usr/bin/env Rscript

library(data.table)
library(matrixStats)

 pi_nuc <- function(vec) {
  alt = sum(vec==1, na.rm=T)
  ref = sum(vec==0,na.rm=T)
  deno = choose(alt+ref,2)
  num = (choose(ref,2) + choose(alt,2))
  pit = (1-num/deno)
  return(pit)
}


site = fread("/plas1/george.sandler/dpgp3/Karl/ZI_allchrom_rawSA_GTs_edited_ord.txt")
site$id = paste(site$V2, site$V3)

targ = fread("/plas1/george.sandler/dpgp3/Karl/male_alt_min5.txt", header = T)
targ$id = paste(targ$V3, targ$V4)
for (i in (1:nrow(targ))){
  df = targ[i,]


  site2 = site[site$id != df$id & site$V2 == df$V3 & site$V3 > (df$V4 -5000)& site$V3 < (df$V4 +5000),]
  site_id = site[site$id == df$id,]

  site3 = rbind(site_id, site2)
  site3[1,194] = "target"

  #site3$id = NULL
  site3$V1 = NULL
  site3$V2 = NULL
  site3$V3 = NULL

  #site4 = t(site3)
  #site4 = dcast(melt(site3, id.vars = "id"), variable ~ id)
  #site4$variable = NULL
  #site4_split <- split(site4, site4$target)

  sit4 = data.frame(t(site3))
  colnames(sit4) <- unlist(sit4[row.names(sit4)=='id',])
  sit4           <- sit4[!row.names(sit4)=='id',]

  #########################################
  #male <- site4[site4$target == 1,]
  #female <- site4[site4$target == 0,]
  
  male <- sit4[as.numeric(as.character(sit4$target)) == 1,]
  female <- sit4[as.numeric(as.character(sit4$target)) == 0,] 
  
  male$target = NULL
  female$target = NULL
 
 
  
  #male <- site4_split$"1"
  #female <- site4_split$"0"
  #########################################
  male <- sapply( male, as.character)
  female <- sapply( female, as.character )
  male <- apply( male,2, as.numeric )
  female <- apply( female,2, as.numeric )

  malmat <- data.matrix(male)
  femalmat <- data.matrix(female)

  mpi = "NULL"
  fpi = "NULL"
  
  mpi= sum(apply(malmat,2,pi_nuc), na.rm=T)/10000
  fpi = sum(apply(femalmat,2,pi_nuc), na.rm=T)/10000
  
  #mpi = mean(rowVars(malmat, na.rm=T))
  #fpi = mean(rowVars(femalmat, na.rm=T))
  #mean(m_s = (rowSums(malmat, na.rm=T)))

  line=(noquote(paste(mpi, fpi)))
  write(line,file="/plas1/george.sandler/dpgp3/Karl/pi_haplotypes_male_alt_min5.txt",append=TRUE)
}


# Jacana SNP Cull for Clinal Loci
source("HZAR_MLE.R")
library(hzar)

#table_cull_hz_hwe <- read.csv("table_cull_hz_hwe_128000.csv",row.names=1)
#locality_info <- read.table("locality_info.txt",header=TRUE)
#str(table_cull_hz_hwe)
#table_to_cull <- cbind(locality_info,table_cull_hz_hwe[as.character(locality_info$Name),])

use.allele="S1_601371118_G"
use.table="table_cull_hz_hwe_128000.csv"
read_cull_table <- function(filename){
  table_cull_hz_hwe <- read.csv(filename,row.names=1)
  locality_info <- read.table("locality_info.txt",header=TRUE)
  table_to_cull <- cbind(locality_info,table_cull_hz_hwe[as.character(locality_info$Name),])
  return(table_to_cull)
}

cull.allele <- function(use.allele,table_cull_hz_hwe) {
  
  use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)

obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                               pObs = table_cull_hz_hwe[[use.allele]],
                               nEff = table_cull_hz_hwe[[use.count]],
                               siteID = table_cull_hz_hwe$locationID)

models=list()
models$free=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "free")
models$fixed=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "fixed")
fitR=lapply(models,hzar.first.fitRequest.old.ML,obsData=obs)
MLE = lapply(fitR,cline.getMLE)

null.AIC=hzar.AICc.hzar.dataGroup(hzar.dataGroup.null(obs))
MLE.AIC=sapply(MLE,hzar.AICc.hzar.cline,nObs=sum(obs$frame$n))

if (any(null.AIC > 2 + MLE.AIC)) return(use.allele)
return(character())
}

if(FALSE){
  #cull.allele(use.allele,read_cull_table(use.table))
  table_to_cull=read_cull_table(use.table)
  allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
  library(doParallel)
  registerDoParallel(2)
  clinal_snps =foreach(use.allele=allele_list,.combine = c)%dopar%{
    cull.allele(use.allele,table_to_cull)
  }
}

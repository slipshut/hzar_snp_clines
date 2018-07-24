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

allele.plot <-function(use.allele,table_cull_hz_hwe, ...) {
  
  use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)
  
  obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                                 pObs = table_cull_hz_hwe[[use.allele]],
                                 nEff = table_cull_hz_hwe[[use.count]],
                                 siteID = table_cull_hz_hwe$locationID)
  hzar.plot.obsData(obs,...)
    
}

allele.cline_plot <-function(use.allele,table_cull_hz_hwe, ...) {
  
  use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)
  
  obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                                 pObs = table_cull_hz_hwe[[use.allele]],
                                 nEff = table_cull_hz_hwe[[use.count]],
                                 siteID = table_cull_hz_hwe$locationID)
  hzar.plot.obsData(obs)
  models=list()
  models$free=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "free")
  models$fixed=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "fixed")
  fitR=lapply(models,hzar.first.fitRequest.old.ML,obsData=obs)
  MLE = lapply(fitR,cline.getMLE)
  
  hzar.plot.cline(MLE$free,add=T,  ...)
  hzar.plot.cline(MLE$fixed,add=T, ...)
  return(MLE)
}
cull.allele <- function(use.allele,table_cull_hz_hwe, plot=F, add = F,MLE.file=NULL,plot.args=list()) {
  
  use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)

obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                               pObs = table_cull_hz_hwe[[use.allele]],
                               nEff = table_cull_hz_hwe[[use.count]],
                               siteID = table_cull_hz_hwe$locationID)

if(!is.null(MLE.file)){
  sel.MLE<-NULL
  if(file.exists(output_file<-paste0(MLE.file,"_",use.allele,".RData"))){
    load(file = output_file)
    if(is.null(sel.MLE) ) return(character(0))
    if(plot) {
      if(sel.MLE$param.all$pMax >0.9 && sel.MLE$param.all$pMin < 0.1){
        if (mean(obs$frame$dist) < (sum(obs$frame$dist * 
                                        obs$frame$obsFreq)/sum(obs$frame$obsFreq)))
          body(sel.MLE$clineFunc) <- bquote(1-.(body(sel.MLE$clineFunc)))
        do.call(hzar.plot.cline,c(list(sel.MLE,add=add),plot.args))
      }
    }
    return(use.allele)
  }
}
models=list()
models$free=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "free")
models$fixed=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "fixed")
fitR=lapply(models,hzar.first.fitRequest.old.ML,obsData=obs)
MLE = lapply(fitR,cline.getMLE)

null.AIC=hzar.AICc.hzar.dataGroup(hzar.dataGroup.null(obs))
MLE.AIC=sapply(MLE,hzar.AICc.hzar.cline,nObs=sum(obs$frame$n))

if (any(null.AIC > 2 + MLE.AIC)) {
  
  if(!is.null(MLE.file)){
    sel.MLE<-MLE[[which.min(MLE.AIC)[1]]]
    save(sel.MLE,file=output_file)
  }
  if(plot) {
    sel.MLE<-MLE[[which.min(MLE.AIC)[1]]]
    if(sel.MLE$param.all$pMax >0.99 && sel.MLE$param.all$pMin < 0.01){
      if (mean(obs$frame$dist) < (sum(obs$frame$dist * 
                                      obs$frame$obsFreq)/sum(obs$frame$obsFreq)))
        body(sel.MLE$clineFunc) <- bquote(1-.(body(sel.MLE$clineFunc)))
      do.call(hzar.plot.cline,c(list(sel.MLE,add=add),plot.args))
    }
  }
  return(use.allele)
}
if(!is.null(MLE.file)){
  sel.MLE<-NULL
  save(sel.MLE,file=output_file)
}
return(character())
}
if(F){
  table_to_cull=read_cull_table(use.table)
  allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
  allele.plot(allele_list[1],table_to_cull, col="transparent")
}

library(doParallel)
registerDoParallel(2)

if(F){
  table_to_cull=read_cull_table(use.table)
  allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
  allele.cline_plot(allele_list[1],table_to_cull)
}

cull.table <- function(use.table,plot=T, add=F, save.MLE=T, plot.args=list()){
  #cull.allele(use.allele,read_cull_table(use.table))
  output_table=sub( x= use.table,pattern="(_[0-9]*\\.csv$)", replacement = "_clinal\\1")
  if(save.MLE)
    output_rdata=file.path("MLE_cache", sub( x= use.table,pattern="(_[0-9]*)\\.csv", replacement = "_clinal\\1"))
  else
    output_rdata=NULL
  if(file.exists(output_table) && (file.mtime(output_table)>file.mtime(use.table))){
    table_to_cull=read_cull_table(output_table)
    allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
    if(!plot) return(allele_list)
    
  } else {
    table_to_cull=read_cull_table(use.table)
    allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
  }
  if(plot){
    if(!add)
      allele.plot(allele_list[1],table_to_cull, col="transparent")
    clinal_snps =foreach(use.allele=allele_list,.combine = c)%do%{
      cull.allele(use.allele,table_to_cull,plot=T,add=T,MLE.file=output_rdata,plot.args=plot.args)
    }
  }else{
    clinal_snps =foreach(use.allele=allele_list,.combine = c)%dopar%{
      cull.allele(use.allele,table_to_cull,plot=F,add=T,MLE.file=output_rdata)
    }
  }
  if(file.exists(output_table) && (file.mtime(output_table)>file.mtime(use.table))){
    return(allele_list)
  }
  use.count=sub(pattern="[ACTG]$",replacement="N",x=clinal_snps)
  col.list=which(colnames(table_to_cull)%in%c(clinal_snps,use.count))
  write.csv(table_to_cull[,c(col.list)],file = output_table)
  return(clinal_snps)
}



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
cull.allele.clinal <- function(use.allele,table_cull_hz_hwe, plot=F, add = F,MLE.file=NULL,
                               plot.args=list(), save.plot=F) {
  
  use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)
  
  obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                                 pObs = table_cull_hz_hwe[[use.allele]],
                                 nEff = table_cull_hz_hwe[[use.count]],
                                 siteID = table_cull_hz_hwe$locationID)
  
  if(!is.null(MLE.file)){
    sel.MLE<-NULL
    if(file.exists(output_file<-paste0(MLE.file,"_",use.allele,".RData"))){
      load(file = output_file)
      if(is.null(sel.MLE) && save.plot) return(list())
      if(is.null(sel.MLE) ) return(character(0))
      #if(plot) {
        if(sel.MLE$param.all$pMax >0.9 && sel.MLE$param.all$pMin < 0.1){
          if (mean(obs$frame$dist) < (sum(obs$frame$dist * 
                                          obs$frame$obsFreq)/sum(obs$frame$obsFreq)))
            body(sel.MLE$clineFunc) <- bquote(1-.(body(sel.MLE$clineFunc)))
          if(save.plot) return(list(sel.MLE$clineFunc))
          if(plot)
            do.call(hzar.plot.cline,c(list(sel.MLE,add=add),plot.args))
          return(use.allele)
        }
      if(save.plot) return(list())
      return(character(0))
      #}
      #return(use.allele)
    }
  }
  stop("Run cull.allele to generate cache files.")
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


cull.hwe.table <- function(use.table,plot=T, add=F, save.MLE=T, plot.args=list(), save.plot = F){
  #cull.allele(use.allele,read_cull_table(use.table))
  output_table=sub( x= use.table,pattern="(_[0-9]*\\.csv$)", replacement = "_clinal\\1")
  if(save.MLE)
    output_rdata=file.path("MLE_cache", sub( x= use.table,pattern="(_[0-9]*)\\.csv", replacement = "_clinal\\1"))
  else
    output_rdata=NULL
  if(file.exists(output_table) && (file.mtime(output_table)>file.mtime(use.table))){
    table_to_cull=read_cull_table(output_table)
    allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
    #if(!plot) return(allele_list)
    
  } else {
    table_to_cull=read_cull_table(use.table)
    allele_list=colnames(table_to_cull)[seq(4,ncol(table_to_cull),by=2)]
  }
  if(save.plot){
    clinal_snps =foreach(use.allele=allele_list,.combine = c)%dopar%{
      cull.allele.clinal(use.allele,table_to_cull,plot=F,add=T,save.plot=T,MLE.file=output_rdata,plot.args=plot.args)
    }
  }else if(plot){
    if(!add)
      allele.plot(allele_list[1],table_to_cull, col="transparent")
    clinal_snps =foreach(use.allele=allele_list,.combine = c)%do%{
      cull.allele.clinal(use.allele,table_to_cull,plot=T,add=T,MLE.file=output_rdata,plot.args=plot.args)
    }
  }else{
    clinal_snps =foreach(use.allele=allele_list,.combine = c)%dopar%{
      cull.allele.clinal(use.allele,table_to_cull,plot=F,add=T,MLE.file=output_rdata)
    }
  }
  if(file.exists(output_table) && (file.mtime(output_table)>file.mtime(use.table))){
    return(clinal_snps)
  }
  #use.count=sub(pattern="[ACTG]$",replacement="N",x=clinal_snps)
  #col.list=which(colnames(table_to_cull)%in%c(clinal_snps,use.count))
  #write.csv(table_to_cull[,c(col.list)],file = output_table)
  return(clinal_snps)
}

# setup plot for clines
setup_cline_plot <- function(...){
  table_cull_hz_hwe<- structure(list(locationID = structure(c(1L, 2L, 3L, 24L, 25L), 
                                                        .Label = c("A","B", "C", "D", "E", "F", "G", "I", "J", "K", "L", "M", "N", "O", 
                                                                                             "P", "Q", "R", "S", "U", "W", "X", "Y", "ZC", "ZD", "ZG"), class = "factor"), 
                                  distance = c(0, 298.1482087, 307.3896796, 672.2232788, 840.8756274
                                  ), 
                                 S1_1070531337_C = c(0.05, 
                                                     0.0625, 0, 1, 0.884615384615385), S1_1070531337_N = c(20L,16L, 2L, 2L, 26L)), .Names = c("locationID", "distance", "S1_1070531337_C", "S1_1070531337_N"), 
                            row.names = c("Palo_Verde", 
                                          "La_Gamba", "Changuinola", "Tocumen", "Aruza_Abajo"), class = "data.frame")
  obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                                 pObs = table_cull_hz_hwe$S1_1070531337_C,
                                 nEff = table_cull_hz_hwe$S1_1070531337_N,
                                 siteID = table_cull_hz_hwe$locationID)
  hzar.plot.obsData(obs,col="transparent",...)
}


# Get cline centers and widths from output
if(!file.exists("Jac_SNP_fixed_hwe_center_width.csv")){
library(foreach)
  
MLE_list <- list.files("MLE_cache",full.names = TRUE, pattern = "table_cull_hz_hwe_clinal_......_S1_([0-9]+)_..RData")
get_snp_cw <- function(MLE_list){
  foreach(x=MLE_list, .combine = rbind) %do% {
  load(x) 
  if(!is.numeric(sel.MLE$param.all$center)) return(data.frame(center=numeric(), width=numeric()))
    if(!(sel.MLE$param.all$pMax >0.9 && sel.MLE$param.all$pMin < 0.1)) 
      return(data.frame(center=numeric(), width=numeric()))
  cat(".")
  data.frame(
  center = sel.MLE$param.all$center,
  width = sel.MLE$param.all$width,
  row.names = sub(replacement="S1_\\1",x=x,
  pattern	= "MLE_cache/table_cull_hz_hwe_clinal_......_S1_([0-9]+)_..RData"))
    }}
 

# Distribution of cline centers
#print(head(snp_cw <- get_snp_cw(sample(MLE_list, 100))))
print(head(snp_cw <- get_snp_cw(MLE_list)))

write.csv(snp_cw, file = "Jac_SNP_fixed_hwe_center_width.csv")
} else{
  snp_cw <- read.csv("Jac_SNP_fixed_hwe_center_width.csv", row.names = TRUE)
}
if(FALSE){
dim(snp_cw) # 8479 this confirms we have 8479 fixed snps

library(lattice)
summary(snp_cw)
histogram(~center, data = snp_cw, nint = 501, xlim=c(400,700), col = "light grey", border = "dark grey", panel = function(...) { panel.histogram(...); 
  panel.abline(v=611,lwd=1,col="purple"); 
  panel.abline(v=c(550.9,672.9),lwd=1,col="purple", lty = "dashed"); 
  panel.abline(v=483.2,lwd=1,col="black");
  panel.abline(v=c(479.2,489.3),lwd=1,col="black", lty="dashed");
  panel.abline(v=490.1, lwd = 1, col = "red")})

densityplot(~center, data = snp_cw, plot.points = FALSE, xlim=c(450,690), lwd = 2, col = "dark grey", panel = function(...) { 
  panel.abline(v=611,lwd=1,col="purple"); 
  panel.abline(v=c(550.9,672.9),lwd=1,col="purple", lty = "dashed"); 
  panel.abline(v=483.2,lwd=1,col="black");
  panel.abline(v=c(479.2,489.3),lwd=1,col="black", lty="dashed"); 
  panel.densityplot(...); 
  panel.abline(v=490.1, lwd = 1, col = "red")})

histogram(~width, data = snp_cw, panel = function(...) { panel.histogram(...); 
  panel.abline(v=288.5,lwd=2,col="purple"); 
  panel.abline(v=c(140.4,485.1),lwd=2,col="purple", lty = "dashed");
  panel.abline(v=32.4,lwd=2,col="black");
  panel.abline(v=c(27.5,41.8),lwd=2,col="black", lty="dashed")})

densityplot(~width, data = snp_cw, plot.points = FALSE, xlim=c(-100,700), lwd = 2, col = "dark grey", panel = function(...) { 
  panel.abline(v=288.5,lwd=2,col="purple"); 
  panel.abline(v=c(140.4,485.1),lwd=2,col="purple", lty = "dashed");
  panel.abline(v=32.4,lwd=2,col="black");
  panel.abline(v=c(27.5,41.8),lwd=2,col="black", lty="dashed");
  panel.densityplot(...); })

densityplot(~log(width), data = snp_cw, plot.points = FALSE, xlim = c(-2,7), lwd = 2, col = "dark grey", panel = function(...) { 
  panel.abline(v=log(288.5),lwd=2,col="purple"); 
  panel.abline(v=log(c(140.4,485.1)),lwd=2,col="purple", lty = "dashed");
  panel.abline(v=log(32.4),lwd=2,col="black");
  panel.abline(v=log(c(27.5,41.8)),lwd=2,col="black", lty="dashed");
  panel.densityplot(...); })

library(moments)
skewness(snp_cw$center) # -17.39349 indicates data is skewed to the left, or negatively skewed
kurtosis(snp_cw$center) # 742.7753 'excess' # leptokurtic - fat-tailed distribution

skewness(snp_cw$width) # 20.5004 indicates data is skewed to the right, or positively skewed
kurtosis(snp_cw$width) # 596.3487

qqnorm(snp_cw$center) # data have more extreme values than would be expected if they truly came from a Normal distribution
  # spike of identical values, heavy tailed

qqnorm(log(snp_cw$width)) # right skewed without log transformation, fat tailed with log transformation

quantile(snp_cw$center, probs=c(0.005,0.995)) # 99% of centers fall within this range
  # 0.5%    99.5% 
#477.0501 542.0634 

skewness(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$center) # 2.365532
  # cut out fat tails to see skewness - cutting out most extreme 1% from both left and right side
kurtosis(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$center) # 12.20143
qqnorm(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$center) # right skewed

quantile(snp_cw$width, probs=c(0.005,0.995)) # 99% of centers fall within this range
  #     0.5%      99.5% 
  # 5.340858 596.137057
skewness(subset(snp_cw, width >= 5.340858 & width <= 596.137057)$width) # 7.595443
  # cut out fat tails to see skewness - cutting out most extreme 1% from both left and right side
kurtosis(subset(snp_cw, width >= 5.340858 & width <= 596.137057)$width) # 102.5963
qqnorm(subset(snp_cw, width >= 5.340858 & width <= 596.137057)$width) # right skewed

skewness(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$width) # 31.22437
# cut out fat tails to see skewness - cutting out most extreme 1% from both left and right side
kurtosis(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$width) # 1191.364
qqnorm(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$width) # right skewed
qqnorm(log(subset(snp_cw, center >= 477.0501 & center <= 542.0634)$width)) # fat tailed
}
## Builiding SNP frequency per locality
#setwd("~/Dropbox/Jacanas/GBS_processing/Formating_Scripts/Diagnostic_loci")
library(foreach)
library(iterators)

for( iter_chunk in 0:37 ){
snp.data <- read.csv(paste0("table_",iter_chunk,".old.csv"))
ind.data <- read.csv("Jacana_251ind_info.csv")

#snp.data[1:10, 1:10]
row.names(snp.data) <- as.character(snp.data$birdA)

#ind.data[1:10, 1:10]
row.names(ind.data) <- as.character(ind.data$IndividualID)

#add locality and site number to snp data (sorted by site), dropping everythign else
snp.data <- cbind(ind.data[,c("Locality","Site_Number")], snp.data[row.names(ind.data),-(1)])

#calculate allele frequency per locality
snp.freq <- foreach(snp.col=colnames(snp.data)[seq(3,ncol(snp.data),2)], .combine=cbind)%do%{
  alleles.col <- sub(pattern = "[ACGT]$", replacement = "N", x=snp.col)
  snp.cnt <- data.frame(snp=tapply(snp.data[[(snp.col)]], snp.data$Locality, sum),
                        tot=tapply(snp.data[[(alleles.col)]], snp.data$Locality, sum))
  snp.cnt$snp <- with(snp.cnt, snp/tot)
  colnames(snp.cnt) <- c(snp.col, alleles.col)
  snp.cnt
}

#summary(which(is.na(snp.freq)))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#       2  3723000  7421000  7414000 11090000 14800000 

#foreach(val=iter(col,by = c("row") ),.combine=rbind ) %:%

#foreach(col=iter(as.matrix(snp.freq),by = c("column") ),.combine=cbind  ) %:% 
#  when(any(is.na(col))) %do% as.data.frame(col[is.na(col),,drop=FALSE])

################################
#error calling combine function:
#  <simpleError in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 8, 17, 5, 15, 11, 1, 16, 2, 13, 12, 3, 4, 10, 6, 7, 18, 14>

#snp.freq[,paste0("S1_185745027_",c("C","N"))]

#print(nextElem(iter(as.matrix(snp.freq),by = c("column"),chunksize=2 )))

snp.freq.cull<-
  foreach(locus.data=iter(as.matrix(snp.freq),by = c("column"),chunksize=2 ),.combine=cbind  ) %:% 
  when(sum(locus.data[,2]>5)>1) %:% when(min(locus.data[locus.data[,2]>5,1])<0.3) %:%
  when(max(locus.data[locus.data[,2]>5,1])>0.7) %:% when(diff(range(locus.data[locus.data[,2]>5,1]))>0.5) %do% 
  as.data.frame(locus.data)

#write out file
#write.csv(snp.freq.cull,file="table_cull_all.csv")

#write out file in pieces
for(chunk in seq(0,ncol(snp.freq.cull)/2-255,by = 250)){
	snp.chunk = snp.freq.cull[,chunk*2+(1:500)]
	write.csv(snp.chunk,file = sprintf("table_cull_hz_%06i.csv",chunk+6000*iter_chunk))
}

chunk = chunk + 250
snp.chunk = snp.freq.cull[,(chunk*2+1):ncol(snp.freq.cull)]
write.csv(snp.chunk,file = sprintf("table_cull_hz_%06i.csv",chunk+6000*iter_chunk))
}

# Jacana SNP Cull for Clinal Loci
source("HZAR_MLE.R")
library(hzar)

table_cull_hz_hwe <- read.csv("table_cull_hz_hwe_128000_test.csv")
str(table_cull_hz_hwe)

use.allele="S1_991670203_T"
use.count=sub(pattern="[ACTG]$",replacement="N",x=use.allele)

obs=hzar.doMolecularData1DPops(table_cull_hz_hwe$distance,
                               pObs = table_cull_hz_hwe[[use.allele]],
                               nEff = table_cull_hz_hwe[[use.count]],
                               siteID = table_cull_hz_hwe$locationID)

models=list()
models$free=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "free")
models$fixed=hzar.makeCline1DFreq(data=obs,tails="none",scaling = "fixed")

MLE = lapply(models,cline.getMLE)

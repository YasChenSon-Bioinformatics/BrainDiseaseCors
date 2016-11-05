source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')

library(annotate)

source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/caprice/R/analyze-GPL570.R')

gpl570_controlled_correlated_subset = list(4523,4522,4358,4218,4136,2821,1917)

all_GDS = list()
j = 1
for (gds in gpl570_controlled_correlated_subset){
  message(paste('GDS', gds, sep=''))
  all_GDS[[j]] = getGEO(GEO=paste('GDS', gds, sep=''), destdir='/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/data')
  j = j + 1
}
all_ESET = list()
j = 1
for (gds in all_GDS){
  all_ESET[[j]] = assign(paste('ESET', all_GDSv[j], sep='_'), GDS2eSet(GDS=gds)) 
  j = j + 1
}
all_DF = list()
j = 1;
for (eset in all_ESET) { 
  all_DF[[j]] = assign(paste('df', all_GDSv[j], sep='_'), as.data.frame(eset))
  j = j + 1
}


all_DF[[4]] = NULL
##########################################################
# Notes:
# Issue with df 4
# SAM
par(mfrow=c(2,3))
par(cex.axis=1.5, cex.lab=1.75, cex.main=1.5, cex.sub=1.5)

bool_to_num = function(x){ if (x==TRUE) 1 else 2}
all_sam = list()
i = 1
for (df in all_DF) {
  # ctrl_count = count(grepl("control", df$disease.state))
  # dz_count = count(!grepl("control", df$disease.state))
  
  dz_ctrl_boolv = grepl('control', df$disease.state)
  y = apply(as.data.frame(dz_ctrl_boolv), 1, bool_to_num)
  
  # y<-c(rep(1,ctrl_count), rep(2,dz_count))
  drops <- c("sample", "age", "gender", "tissue", "genotype.variation", "development.stage", "agent", "other", "cell.type", "disease.state", "description", "individual")
  df = df[, !(names(df) %in% drops)]
  
  df_t = t(df)
  
  samfit <- SAM(df_t, y, resp.type="Two class unpaired", nperms=1)
  message(i, ' -----------------------------------------------')
  message(y)
  plot(samfit)
  all_sam[[i]] = samfit
  i = i + 1
}


all_genes = list()
# Put up and down regulated genes in a list of lists
i = 1
for (sam in all_sam) {
  all_genes[[i]] = list(all_sam[[i]]$siggenes.table$genes.up, all_sam[[i]]$siggenes.table$genes.lo)
  i = i + 1
}


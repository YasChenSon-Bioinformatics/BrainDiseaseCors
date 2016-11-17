source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')

library(annotate)

# source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/caprice/R/analyze-GPL570.R')

gpl570_controlled_correlated_subset = list(4523,4522,4358,4218,4136,2821,1917)
gds_short = gpl570_controlled_correlated_subset

all_GDS = list()
j = 1
for (gds in gpl570_controlled_correlated_subset){
  message(paste('GDS', gds, sep=''))
  all_GDS[[j]] = getGEO(GEO=paste('GDS', gds_short[j], sep=''), destdir='/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/data')
  j = j + 1
}
all_ESET = list()
j = 1
for (gds in all_GDS){
  all_ESET[[j]] = assign(paste('ESET', gds_short[j], sep='_'), GDS2eSet(GDS=gds)) 
  j = j + 1
}
all_DF = list()
j = 1;
for (eset in all_ESET) { 
  all_DF[[j]] = assign(paste('df', gds_short[j], sep='_'), as.data.frame(eset))
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
  
  samfit <- SAM(df_t, y, resp.type="Two class unpaired", nperms=100)
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


library(gdata)
# Check results for 4358
supp_4358 = read.xls('~/Desktop/Columbia/Bioinformatics/project/paper DEGs data/4358_supp.xls')
genes_4358 = c(37872,23322,19301,13727,13887,11756,48990,18455,26742,31259,11879,18598,18106,27227,13019,39660,27829,12271,19170,11719,13316,9038,22273,34673,19376,11050,32499,30316,19183,10353,23753,10609,23431,15779,13699,14268,25817,10354,22595,28271,24608,14254,13518,11718,19253,31753,22044,19179,11642,34092,31156,18228,23232,32280,36355,18397,45608,23101,40832,36436,31097,11895,48255,35731,15538,21116,13421,22355,12890,31907,18382,37514,36483,10764,2622,11582,29148,30803,21720,18553,11860,21220,20761,18544,3533,18555,18145,14545,23759,44893,2607,41211,994,15580,10435,24464,35757,23974,14931,26783,33783,37407,33961,12494,11827,34218,51484,28869,27025,15914,10336,32073,11535,51166,21294,35877,14995,19061,24416,52521,35463,54632,51316,24304,11806,45689,28494,18377,10365,5810,36351,53596,11885,12312,20855,39255,2418,36748,13302,36624,39688,18163,26211,21384,35187,27416,21919,22805,27106,12682,37493,38705,14354,8843,38449,33822,14297,35211,17856,25856,17442,40140,12956,18365,54309,14377,22151,14269,21592,35189,14195,54527,28790,20856,7552,33946,23503,37326,28804,34431,33922,11488,18901,11886,11515,1272,44671,13171,21800,24105,13916,39990,39072,54623,37415,13412,28637,12597,23359,11247,21743,53288,18800,48011,19826,21254,53038,15820,28174,21510,19905,45106,33126,31230,35258,54651,41630,53843,52664,49229)
overlap = supp_4358[which(supp_4358$GeneName %in% genes_4358),]



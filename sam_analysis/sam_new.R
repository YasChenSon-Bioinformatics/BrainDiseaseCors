source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')
library('limma')


source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis/functions.R')
setwd('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis')

datasets_num = list(4523, 4522, 4358, 4218, 4136, 2821, 1917, 5204, 4135, 2795, 1962)
datasets = download_GDS(datasets_num)



# COMPUTE SAM ##########################################################################
bool_to_num = function(x){ if (x==TRUE) 1 else 2}
all_sam = list()
i = 1
for (df in datasets) {
  message(datasets_num[[i]], ' -----------------------------------------------')
  
  if (i < 4) { i = i + 1; next; } # uncomment to fast forward
  if (i > 4) { i = i + 1; next; } # uncomment to fast forward
  # if (i == 4) { i = i + 1; next } # error in 4th dataset...
  
  # Possible levels used for control are:
  # 
  # "control"
  # "normal" (2795)
  # "young" (aging study 5204)
  # "Braak stage I-II" (alz 4135)
  # "non-tumor" (scf/angiogenesis 1962)
  
  control_strings = c("control", "normal", "Braak stage I-II", "non-tumor")
  dz_ctrl_boolv = grepl(paste(control_strings, collapse = "|"), df$disease.state)
  
  # Aging study has no disease.state column - use "age"
  # Note - grepl doesn't seem to like special characters
  if (datasets_num[[i]] == 5204) {
    age_control_strings = c("young", "middle")
    dz_ctrl_boolv = grepl(paste(age_control_strings, collapse = "|"), df$age)
  }
  
  y = apply(as.data.frame(dz_ctrl_boolv), 1, bool_to_num)
  
  drops <- c("sample", "age", "gender", "tissue", "genotype/variation", "genotype.variation", "development.stage", "agent", "other", "cell.type", "disease.state", "description", "individual")
  df = df[, !(names(df) %in% drops)]
  df_t = t(df)
  
  samfit <- SAM(df_t, y, resp.type="Two class unpaired", nperms=2, fdr.output = 0.01, geneid = rownames(df_t))
  all_sam[[i]] = samfit
  i = i + 1
}
#########################################################################################




# PLOT SAM ##############################################################################
par(mfrow=c(3,3))
par(cex.axis=1.5, cex.lab=1.75, cex.main=1.5, cex.sub=1.5)
for (sam in all_sam) {
  plot(sam)
}


# Collect DEGs ##############################################################################
i = 1
deg_up = list()
deg_lo = list()
for (sam in all_sam) {
  message(i, typeof(all_sam[[i]]$siggenes.table$genes.up), typeof(all_sam[[i]]$siggenes.table$genes.lo))
  if (i == 4) { i = i + 1; next } # error in 4th dataset...
  
  # Take the Gene Names (e.g. 6473)
  if (typeof(all_sam[[i]]$siggenes.table$genes.up) == "matrix") {
    deg_up[[datasets_num[[i]]]] = all_sam[[i]]$siggenes.table$genes.up[,2] 
  } else if (typeof(all_sam[[i]]$siggenes.table$genes.up) == "character") {
    deg_up[[datasets_num[[i]]]] = all_sam[[i]]$siggenes.table$genes.up[,2]
  }
  
  if (typeof(all_sam[[i]]$siggenes.table$genes.lo) == "matrix") {
    deg_lo[[datasets_num[[i]]]] = all_sam[[i]]$siggenes.table$genes.lo[,2] 
  } else if (typeof(all_sam[[i]]$siggenes.table$genes.lo) == "character") {
    deg_lo[[datasets_num[[i]]]] = all_sam[[i]]$siggenes.table$genes.lo[,2] 
  }
  
  i = i + 1
}

# Save CSVs
for (i in 1:length(deg_up)) {
  write.table(deg_up[[i]], paste(datasets_num[[i]], '_genes_up.csv', sep=''), sep=',')
  write.table(deg_lo[[i]], paste(datasets_num[[i]], '_genes_lo.csv', sep=''), sep=',')         
}




# Count DEGs per dataset
library(dplyr)
setwd('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/caprice/RESULT')
sam_results = read.csv('11GDS_sam_results.csv')
counts = sam_results[!duplicated(sam_results$gds),] %>% select(gds, count)

# gds count
# <int> <int>
#   1   4522     4
# 2   4358   850
# 3   4218  5372
# 4   2821    34
# 5   5204  8663
# 6   4135 28186
# 7   1962 27691
# 8   4523     2
# 9   4136    33
# 10  1917    19
# 11  2795     3

boxplot(count ~ gds, counts)
boxplot(log(count) ~ gds, counts)
boxplot(count ~ gds, counts, log='y', ylab='log DEG count')


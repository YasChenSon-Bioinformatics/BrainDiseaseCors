source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')
library('limma')


source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis/functions.R')
setwd('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis')

datasets_num = list(1917, 1962, 2795, 2821, 4522, 4523, 4135, 4136, 4218, 4358, 5204)
datasets = download_GDS(datasets_num)



# COMPUTE SAM (two class dz/ctrl) ##########################################################################
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




# PLOT SAM ##############################################################################
plot_sam <- function(sam_list, gds_num) {
  par(mfrow=c(3,3))
  par(cex.axis=1.5, cex.lab=1.75, cex.main=1.5, cex.sub=1.5)
  i = 1
  for (sam in sam_list) {
    plot(sam)
    title(main = paste0(gds_num[[i]]))
    i = i + 1
  }
}


# Collect DEGs ##############################################################################
collect_degs <- function(all_sam) {
  
  i = 1
  deg_up = list()
  deg_lo = list()
  
  for (sam in all_sam) {
    message(i, typeof(all_sam[[i]]$siggenes.table$genes.up), typeof(all_sam[[i]]$siggenes.table$genes.lo))
    
    # Take the Gene Names (e.g. 6473)
    deg_up[[ as.character(datasets_num[[i]]) ]] = all_sam[[i]]$siggenes.table$genes.up[,2]
    deg_lo[[ as.character(datasets_num[[i]]) ]] = all_sam[[i]]$siggenes.table$genes.lo[,2] 
    
    i = i + 1
  }
  
  deg = list()
  deg[['lo']] = deg_lo
  deg[['up']] = deg_up
  return (deg)
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



# Multiclass ####################################################################
#
# 5204 young - middle - normal - extremely aged
# 4523 control - schizophrenia
# 4522 control - schizophrenia
# 4358 control - HIV - HIV + HAD - HIV + HAD + HIVE
# 4218 healthy control - ms early stage active inflammation - ms after demylenation active inflammation - ms after inflammation late stage
# 4136 control - incipient - moderate - severe stage
# 4135 Braak stage I/II - Braak stage III/IV - Braak stage V/VI
# 2821 control - parkinson's disease
# 2795 normal - neurofibrillary tangle
# 1962 non-tumor - astrocytomas - oligodendrogliomas - glioblastomas
# 1917 control - schizophrenia
multiclass_df_num = c(1962, 4135, 4136, 4358, 5204) # 4218 - only one sample with "MS-early stage-active inflammation" - breaks SAM
multiclass_df = datasets[which(datasets_num %in% multiclass_df_num)]

# Copy and paste these instead of looking up columns by typing $ (takes forever)
# multiclass_df[[2]]$disease.state
# datasets[[5]]$disease.state

i = 0
multiclass_samfit = list()
for (df in multiclass_df) {
  i = i + 1
  message(multiclass_df_num[[i]], ' -----------------------------------------------')
  # Compute numeric class labels for samples
  # based on disease level
  y <- as.factor(df$disease.state)
  
  # browser()
  
  # Aging study has no disease.state column - use "age"
  # Note - grepl doesn't seem to like special characters
  if (multiclass_df_num[[i]] == 5204) {
    y <- as.factor(df$age)
  }
  
  levels(y) <- 1:length(levels(y))
  y <- as.numeric(y)
  
  drops <- c("sample", "age", "gender", "tissue", "genotype/variation", "genotype.variation", "development.stage", "agent", "other", "cell.type", "disease.state", "description", "individual")
  df = df[, !(names(df) %in% drops)]
  df_t = t(df)
  
  samfit <- SAM(df_t, y, resp.type="Multiclass", nperms=2, fdr.output = 0.01, geneid = rownames(df_t))
  multiclass_samfit[[i]] = samfit
}

plot_sam(multiclass_samfit, multiclass_df_num)
degs = collect_degs(multiclass_samfit)
degs[[1]]





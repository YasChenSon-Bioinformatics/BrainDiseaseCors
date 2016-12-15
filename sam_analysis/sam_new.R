source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')
library('limma')


source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis/functions.R')
setwd('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis')

datasets_num = list(1917, 1962, 2795, 2821, 4522, 4523, 4135, 4136, 4218, 4358, 5204)
gds_esets_dfs = download_GDS(datasets_num)
gds = gds_esets_dfs$all_GDS
esets = gds_esets_dfs$all_ESET
dfs = gds_esets_dfs$all_DF


# COMPUTE SAM (two class dz/ctrl) ##########################################################################
# @param random - Randomly permute the sample labels to generate random list of DEGs from SAM for comparison to actual SAM
# Note - too many N/A for 4218!
bool_to_num = function(x){ if (x==TRUE) 1 else 2}

do_sam = function(gds, dfs, random) {
  
  all_sam = list()
  
  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    expCond <- gds[[i]]@dataTable@columns
    
    message(datasets_num[[i]], ' -----------------------------------------------')
    
    # if (i < 8) { next; } # uncomment to fast forward
    # if (i == 4) {
    
    # Possible levels used for control are:
    # 
    # "control"
    # "normal" (2795)
    # "young" (aging study 5204)
    # "Braak stage I-II" (alz 4135)
    # "non-tumor" (scf/angiogenesis 1962)
    
    control_strings = c("control", "normal", "Braak_stage_I-II", "non-tumor")
    dz_ctrl_boolv = grepl(paste(control_strings, collapse = "|"), expCond$disease.state)
    
    # Aging study has no disease.state column - use "age"
    # Note - grepl doesn't seem to like special characters
    if (datasets_num[[i]] == 5204) {
      age_control_strings = c("young", "middle")
      dz_ctrl_boolv = grepl(paste(age_control_strings, collapse = "|"), expCond$age)
    }
    
    y = apply(as.data.frame(dz_ctrl_boolv), 1, bool_to_num)
    
    # Use random boolv
    if (random == TRUE) {
      y <- sample(1:2, length(y), replace=T)
    }
    
    # FIXME: change nperm from 2 (debug) to 1000 (production)
    samfit <- SAM(df, y, resp.type="Two class unpaired", nperms=1000, fdr.output = 0.1,
                  geneid = rownames(df) )
    all_sam[[i]] = samfit
  }

  return(all_sam)

}

sam_real = do_sam(datasets)
sam_random = do_sam(gds, dfs, random = TRUE)



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
  
  samfit <- SAM(df_t, y, resp.type="Multiclass", nperms=2, fdr.output = 0.01, genename = rownames(df_t))
  multiclass_samfit[[i]] = samfit
}

plot_sam(multiclass_samfit, multiclass_df_num)
degs = collect_degs(multiclass_samfit)
degs[[1]]



# DEG network ################################################
#
#
library(ggnetwork)
library(network)
library(sna)

compute_network <- function(all_sam) {
  deg <- list()
  
  type_up <- 'up'
  type_lo <- 'lo'
  
  library(dplyr)
  uped <- rbindlist(
    sapply(all_sam, function(sam) sam$siggenes.table$genes.up %>% as.data.frame )[
      which(sapply(all_sam, function(sam) ! is.null(sam$siggenes.table$genes.up)))
      ])
  uped$gds = rep(unlist(datasets_num), sapply(all_sam, function(x) x$siggenes.table$ngenes.up))
  uped$type = type_up
  
  loed <- rbindlist(
    sapply(all_sam, function(sam) sam$siggenes.table$genes.lo  %>% as.data.frame )[
      which(sapply(all_sam, function(sam) ! is.null(sam$siggenes.table$genes.lo)))
      ])
  loed$gds = rep(unlist(datasets_num), sapply(all_sam, function(x) x$siggenes.table$ngenes.lo))
  loed$type = type_lo
  
  degdf <- bind_rows(uped, loed)
  colnames(degdf) <- c("geneid", 'probe', 'score', 'numerator',
                       'denominator', 'fold-change', 'qval', 'gds', 'type')
  write.csv(degdf, file = "caprice/RESULT/11GDS_sam_results.csv", quote=FALSE, row.names=FALSE)
  
  degdf %>% group_by(probe, type) %>% mutate( n_dataset = n() ) %>% filter( n_dataset > 1 ) %>% summarize( gds = collapse() )
  
  degdf %>% group_by(gds, type) %>% summarize(n = n() )
  
  k <- 1
  i <- 2
  
  deg_matrix <- list()
  
  deg_matrix <- matrix( NA, ncol = length(datasets), nrow = length(datasets) )
  
  rownames(deg_matrix) <- paste0('GDS', datasets_num) 
  colnames(deg_matrix) <- paste0('GDS', datasets_num)
  
  for( i in seq_along(datasets) ){
    for( k in i:length(datasets) ){
      
      gds_i <- datasets_num[[i]]
      gds_k <- datasets_num[[k]]
      
      deg_i_up <- (degdf %>% filter( gds == gds_i & type == type_up ))$probe
      deg_i_lo <- (degdf %>% filter( gds == gds_i & type == type_lo ))$probe
      deg_k_up <- (degdf %>% filter( gds == gds_k & type == type_up ))$probe
      deg_k_lo <- (degdf %>% filter( gds == gds_k & type == type_lo ))$probe
      
      if( i == k )
        next
      deg_matrix[i,k] <- length(intersect(deg_i_up, deg_k_up))
      deg_matrix[k,i] <- length(intersect(deg_i_lo, deg_k_lo))
    }
  }
  
  #         GDS4523 GDS4522 GDS4358 GDS4218 GDS4136 GDS2821 GDS1917 GDS5204 GDS4135 GDS2795 GDS1962
  # GDS4523      NA       0       0       0       0       0       0       0       0       0       0
  # GDS4522       0      NA       0       0       0       0       0       0       0       0       0
  # GDS4358       0       0      NA       0       0       0       0       0       0       0       0
  # GDS4218       0       1       0      NA       0       0       0       0       0       0       0
  # GDS4136       0       0       0       0      NA       0       0       0       0       0       0
  # GDS2821       0       0       0       1       0      NA       0       0       0       0       0
  # GDS1917       0       0       0       1       0       0      NA       0       0       0       0
  # GDS5204       0       1       0     564       0       2       2      NA       0       0       0
  # GDS4135       0       0       0       0       0       0       0       0      NA       0       0
  # GDS2795       0       0       0       0       0       0       0       0       0      NA       0
  # GDS1962       0       1       0    1054       0       5       2    4007       0       0      NA
  
  
  my_deg <- sapply(all_sam, function(sam) sam$siggenes.table$ngenes.up + sam$siggenes.table$ngenes.lo )
  
  dnet <- network(deg_matrix, #vertex.attr = my_deg, vertex.attrnames = rep('thisDEG',length(my_deg)),
                  directed = FALSE, bipartite = FALSE)
  dnet <- set.edge.value(dnet, 'n_DEG_lo',   deg_matrix)
  dnet <- set.edge.value(dnet, 'n_DEG_up', t(deg_matrix))
  as.matrix(dnet,attrname='n_DEG_lo')
  as.matrix(dnet,attrname='n_DEG_up')
  
  set.seed(1234)
  png("caprice/RESULT/sam-deg-network.png", width = 960, height=960)
  ggplot(dnet, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50", aes(size=n_DEG_up+n_DEG_lo) ) +
    geom_edgetext_repel(aes(label = n_DEG_up), color = "black", size = 7,
                        fill = "#E6AAAA", box.padding = unit(1.5, "lines") ) +
    geom_edgetext_repel(aes(label = n_DEG_lo), color = "black", size = 7,
                        fill = "skyblue", box.padding = unit(.5, "lines")) +
    geom_nodes(color = "white", size = 17) +
    geom_nodetext(aes(label = paste0("GDS\n",substr(vertex.names,4,100)) ),
                  fontface = "bold", color="black", size=8, lineheight=.8) +
    theme_blank()
  dev.off()
}


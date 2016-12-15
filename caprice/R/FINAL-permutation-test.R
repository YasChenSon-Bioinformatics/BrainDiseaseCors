source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
Ml <- extractMatrixFromEset(ESETl)

non_permutated <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=500, removeGenderEffect=TRUE,method='fdr', permute=NULL)
true_deg <- build_deg_matrix(non_permutated)
for( k in seq_along(gdsv) ){
    peal[[k]] <- do_PathwayEnrichmentAnalysis(gene2pathwaydf, non_permutated, gdsv, gds = gdsv[k])
}
pea_matrix <- do.call(cbind, lapply( peal, function(x) x$pval ))
colnames(pea_matrix) <- gdsv
rownames(pea_matrix) <- peal[[1]]$pathway
pea_bitmatrix <- pea_matrix < arbitrary_p_threshold 
true_common_pea_matrix <- t(pea_bitmatrix) %*% pea_bitmatrix

gene2pathwaydf <-
    #read.delim("caprice/MAP/UniProt2Reactome_All_Levels.tsv",  # too minute
    read.delim("caprice/MAP/UniProt2Reactome.tsv",
               sep='\t', header = TRUE, stringsAsFactors=FALSE) %>%
    mutate( pathway = gsub(' ','_',pathway) )


gdsv <- sapply(GDSl, function(g) g@header$dataset_id[1])
g <- GDSl[[1]]
out <- list()
deg_mat <- list()
peal <- list() # pea list
common_pea_matrices <- list()
arbitrary_p_threshold  <- .05
for( i in 1:100 ){
    message(i, "\n")
    set.seed(i)
    out <- suppressMessages(
        applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=500, removeGenderEffect=TRUE, permute=TRUE)
    )
    deg_mat[[i]] <- build_deg_matrix(out)
    for( k in seq_along(gdsv) ){
        peal[[k]] <- do_PathwayEnrichmentAnalysis(gene2pathwaydf, out, gdsv, gds = gdsv[k])
    }
    pea_matrix <- do.call(cbind, lapply( peal, function(x) x$pval ))
    colnames(pea_matrix) <- gdsv
    rownames(pea_matrix) <- peal[[1]]$pathway
    pea_bitmatrix <- pea_matrix < arbitrary_p_threshold 
    common_pathway_matrix <- t(pea_bitmatrix) %*% pea_bitmatrix
    common_pea_matrices[[i]] <- common_pathway_matrix
}


deg_mat[[1]]


sapply(deg_mat, function(m) m['GDS1962', 'GDS5204'] ) > 

sapply(deg_mat, function(m) m['GDS5204', 'GDS1962'] )
sapply(deg_mat, function(m) m['GDS4522', 'GDS4523'] )
sapply(deg_mat, function(m) m['GDS4523', 'GDS4522'] )

save(deg_mat, file='MISC/deg_mat.Rdata')
deg_list <- lapply(deg_mat, function(m) (m - true_deg) >= 0)

pval_matrix <- matrix(data=0, nrow = length(GDSl), ncol = length(GDSl) )
for( i in 1:1000) {
    pval_matrix <- pval_matrix + deg_list[[i]]
}
pval_matrix <- pval_matrix * (1/1000)
#         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# GDS5204      NA   1.000   1.000   0.267   1.000   0.746   1.000   1.000   0.622   0.638   1.000
# GDS4522   0.558      NA   0.002   0.104   0.687   0.785   0.391   0.010   0.528   0.149   0.581
# GDS4523   0.565   0.004      NA   0.090   0.205   0.779   0.120   0.013   0.193   0.367   1.000
# GDS4358   1.000   1.000   1.000      NA   1.000   0.039   0.122   1.000   0.045   0.105   0.246
# GDS4218   1.000   0.009   1.000   1.000      NA   1.000   1.000   0.735   0.694   0.692   0.688
# GDS4136   0.000   0.387   1.000   1.000   1.000      NA   0.029   0.554   0.145   0.299   0.239
# GDS4135   1.000   1.000   0.539   1.000   1.000   1.000      NA   0.284   0.006   0.261   0.017
# GDS2821   0.010   0.139   0.278   1.000   0.159   0.169   1.000      NA   0.699   0.228   0.739
# GDS2795   0.037   1.000   1.000   1.000   1.000   1.000   1.000   0.581      NA   0.230   0.674
# GDS1962   0.000   1.000   0.278   1.000   1.000   0.048   0.463   1.000   0.521      NA   0.343
# GDS1917   0.606   0.528   0.663   1.000   0.511   0.520   1.000   0.612   0.655   1.000      NA


path_matrix <- matrix(data=0, nrow = length(GDSl), ncol = length(GDSl) )
for( i in 1:100) { # pathway computation is really heavy ... just 100 times
    bit_matrix <- (common_pea_matrices[[i]] - true_common_pea_matrix) >= 0
    path_matrix <- path_matrix + bit_matrix
}
path_matrix <- path_matrix * (1/100)
#         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# GDS5204    0.09    1.00    1.00    0.09    0.06    0.31    0.57    0.15    0.42    0.00    0.18
# GDS4522    1.00    1.00    0.27    0.44    0.66    0.09    1.00    1.00    1.00    1.00    0.54
# GDS4523    1.00    0.27    0.62    1.00    1.00    0.58    0.51    0.14    0.06    1.00    0.13
# GDS4358    0.09    0.44    1.00    0.56    0.31    1.00    0.18    0.28    0.17    0.45    0.06
# GDS4218    0.06    0.66    1.00    0.31    0.61    0.62    1.00    0.61    1.00    0.00    0.57
# GDS4136    0.31    0.09    0.58    1.00    0.62    0.29    1.00    0.51    0.61    0.18    1.00
# GDS4135    0.57    1.00    0.51    0.18    1.00    1.00    0.78    0.66    0.21    0.30    0.55
# GDS2821    0.15    1.00    0.14    0.28    0.61    0.51    0.66    0.26    0.20    0.08    0.71
# GDS2795    0.42    1.00    0.06    0.17    1.00    0.61    0.21    0.20    0.02    0.00    0.09
# GDS1962    0.00    1.00    1.00    0.45    0.00    0.18    0.30    0.08    0.00    0.02    0.58
# GDS1917    0.18    0.54    0.13    0.06    0.57    1.00    0.55    0.71    0.09    0.58    0.05

# remarkable result between 1962 and 5204
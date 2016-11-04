source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

# search Differentially Expressed Genes for each datasets 

GDSl <- download_GDSs(GDSnumberv = c(4523, 4522, 4358, 4218, 4136, 2821, 1917), skipv = c('not-GPL570'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl)
Ml <- extractMatrixFromEset(ESETl)

plotDegCorMatrix <- function( Ml, GDSl, n = 400 ){
  topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene = n, method = 'fdr') # CAUTION: currently our design matrices are just 'disease.state'. no age or gender
  DEG_matrix <- sapply(topped, function(x) rownames(x$table))
  length(       as.vector(DEG_matrix) ) # 2800
  length(unique(as.vector(DEG_matrix))) # 2720
  
  # FIXME: need to deal with GPL570 and hgu133a difference
  
  commonDEGv <- names(table(DEG_matrix)[table(DEG_matrix) > 1])
  GDSnamev   <- sapply(GDSl, function(gds){gds@header$dataset_id[1]})
  
  i <- 2
  disease_signature_matrix <- matrix(NA, nrow = length(commonDEGv), ncol = length(ESETl), dimnames =list(commonDEGv, GDSnamev) )
  rownames(disease_signature_matrix)
  for ( i in seq_along(ESETl) ){
    this <- exprs(ESETl[[i]])[commonDEGv, ]
    control_samplev <- grepl('control',GDSl[[i]]@dataTable@columns$disease.state)
    disease_signature_matrix[ , i ] <-  rowMeans(this[ ,   control_samplev], na.rm=TRUE) -
      rowMeans(this[ , ! control_samplev], na.rm=TRUE)
    # NB: GDS4522, GDS4218, GDS2821 is not on GPL570 so na.rm=FALSE is necessary
  }
  # if all columns of a particular proble (e.g. 202018_s_at in GDS) are NA, rowMeans returns NaN
  #   i <- 2 ; exprs(ESETl[[i]])[rownames(exprs(ESETl[[i]]))=='202018_s_at',]
  #   i <- 4 ; exprs(ESETl[[i]])[rownames(exprs(ESETl[[i]]))=='1559746_a_at',]
  #
  # cor cannot handle NaN but can handle NA, thus convert
  disease_signature_matrix[is.nan(disease_signature_matrix)] <- NA
  
  #disease_signature_matrix <- as.matrix(na.omit(as.data.frame(disease_signature_matrix))) # remove rows with na 

  source('~/BrainDiseaseCors/caprice/R/pairs-extension.R')
  pairs(disease_signature_matrix, lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,
        main = paste0("Matrix of Correlation Coefficients\n", "number of DEGs for each dataset =", n) )
}

pdf('~/BrainDiseaseCors/caprice/EDA/nDeg-affects-cors.1stEd.pdf', width=15,height=12)
for ( n in c(200, 300, 400, 500, 600, 1000, 1500, 3000, 6000) ){
  plotDegCorMatrix(Ml, GDSl, n = n)
}
dev.off()
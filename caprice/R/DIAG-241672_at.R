setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(c(1962, 5204), skipv='');
ESETl <- convertGDS2ESET(GDSl)
Ml <- extractMatrixFromEset(ESETl)


out <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene = 100)
out

intersect(rownames(out[[1]]$table), rownames(out[[2]]$table))

p <- '241672_at'

par(mfrow=c(1,2))
plot(Ml[[1]][p, ], col=as.factor(GDSl[[1]]@dataTable@columns$disease.state=='non-tumor'), main=paste0(p, ' in GDS1962 (Tumor Study)\n Red: non-tumor'))
plot(Ml[[2]][p, ], col=as.factor(GDSl[[2]]@dataTable@columns$age%in%c('young (<40yr)', 'middle aged (40-70yr)')), main=paste0(p, ' in GDS5204 (Aging study)\n Red: non-elderly(<70yr)'))

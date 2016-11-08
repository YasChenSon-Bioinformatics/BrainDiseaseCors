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
plot(Ml[[1]][p, ], lwd=3, col=as.factor(GDSl[[1]]@dataTable@columns$disease.state=='non-tumor'),
     pch=as.numeric(GDSl[[1]]@dataTable@columns$disease.state=='non-tumor')+2,
     main=paste0(p, ' in GDS1962 (Tumor Study)\n Red-Plus: non-tumor'), ylab="Expression Values")
plot(Ml[[2]][p, ], col=as.factor(GDSl[[2]]@dataTable@columns$age%in%c('young (<40yr)', 'middle aged (40-70yr)')), 
     ch=as.numeric(GDSl[[2]]@dataTable@columns$age%in%c('young (<40yr)', 'middle aged (40-70yr)'))+2,
     main=paste0(p, ' in GDS5204 (Aging study)\n Red-Plus: non-elderly(<70yr)'), ylab="Expression Values", lwd=3)

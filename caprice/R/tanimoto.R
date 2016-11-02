source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')

setGlobalConstantList()
loadLibraries()

threeGDSl <- download_GDSs(GDSnumberv = c(4602, 4491, 4589), skipv=c())
ESETl <- convertGDS2ESET(threeGDSl)
Ml <- extractMatrixFromEset(ESETl)


plotMA( ESETl[[1]], array=1, xlab="Average log-expression (A-values)", ylab="Expression log-ratio (this sample vs others, M-values)")
plotMA( ESETl[[2]], array=1, xlab="Average log-expression (A-values)", ylab="Expression log-ratio (this sample vs others, M-values)")
plotMA( ESETl[[3]], array=1, xlab="Average log-expression (A-values)", ylab="Expression log-ratio (this sample vs others, M-values)")

count <- 0
png(paste0('BrainDiseaseCors/caprice/EDA/MAplots-ofThreeGDS.1stEd.', count / 16 ,'.png'), width = 1000, height=1000)
for ( i in seq_along(ESETl) ){
  #for ( i in 1:2 ){
  for ( k in seq_along(colnames(ESETl[[i]]@assayData$exprs))  ){
    if (count %% 16 == 0 ){
      dev.off()
      png(paste0('BrainDiseaseCors/caprice/EDA/MAplots-ofThreeGDS.1stEd.', count / 16 ,'.png'), width = 1000, height=1000)
      par(mfrow=c(4,4))
    }
    count <- count + 1
    p <- plotMA( ESETl[[i]], array=k,
                 xlab="Average log-expression (A-values)",
                 ylab="Expression log-ratio (this sample vs others, M-values)",
                 main = paste0('DATASET ', i, ' ', threeGDSl[[i]]@header$dataset_id[1],
                               '\n (Sample ',k,') ', colnames(ESETl[[i]]@assayData$exprs)[k]))    
    abline(h = 0, col='blue',   lwd=2, lty=2)
    abline(v = 1, col='red',    lwd=2, lty=2)        
    abline(v = 2, col='blue',   lwd=2, lty=2)
    abline(v = 3, col='green',  lwd=2, lty=2)        
    abline(v = 4, col='orange', lwd=2, lty=2)        
    print(p)
  }
  message(i)
}
dev.off()

# convert *.png MAplots-ofGDSsOnGPL570.1stEd.pdf
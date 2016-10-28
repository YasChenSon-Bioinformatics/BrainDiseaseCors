
source('/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/R/analyze-GPL570.R')

setGlobalConstantList(rootDir = '/home/PCUser/BrainDiseaseCors')
loadLibraries()
allGDSl <- donwload_GDSs(skipv=c('not-GPL570'))

#GDS_on_GPL570v <- sapply( all_GDS, function(gds) {gds@header$dataset_id[1]} )

ESETl <- convertGDS2ESET(allGDSl)

#pdf('/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/EDA/MAplots-ofGDSsOnGPL570_.1stEd.pdf', width=10, height=10)

count <- 0
for ( i in seq_along(ESETl) ){
#for ( i in 1:2 ){
    for ( k in seq_along(colnames(ESETl[[i]]@assayData$exprs))  ){
        if (count %% 16 == 0 ){
            dev.off()
            png(paste0('/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/',
                       'BrainDiseaseCors/caprice/EDA/MAplots-ofGDSsOnGPL570.1stEd.', count / 16 ,'.png'),width = 1000, height=1000)
            par(mfrow=c(4,4))
        }
        count <- count + 1
        p <- plotMA( ESETl[[i]], array=k,
                xlab="Average log-expression (A-values)",
                ylab="Expression log-ratio (this sample vs others, M-values)",
                main = paste0('DATASET ', i, ' ', allGDSl[[i]]@header$dataset_id[1],
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
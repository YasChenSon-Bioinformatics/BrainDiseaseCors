example <- function(){
    library(RankProd)
    
    data(arab) # normalized by RMA (i.e. log2 scale)
    head(arab)
    glimpse(arab)
    boxplot(arab)
    arab.cl # class label
    arab.origin # lab Number
    
    arab.gnames # gene names (AffyID)
    
    
    # RankProducts() is simpler.
    # RP.advance() for multiple origins.
    
    rpout <-
        RankProducts(
            data = arab[ , which(arab.origin==1)],
            cl = arab.cl[which(arab.origin==1)],
            gene.names = arab.gnames,
            logged = TRUE, # already log2 scale
            na.rm=FALSE, plot=FALSE, rand=123)
}



source('/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/R/analyze-GPL570.R')

setGlobalConstantList()
loadLibraries()
allGDSl <- download_GDSs(skipv=c('not-GPL570', 'no-disease.state', 'no-control',
                                 'non-binary', 'found-NA'))

#GDS_on_GPL570v <- sapply( all_GDS, function(gds) {gds@header$dataset_id[1]} )

ESETl <- convertGDS2ESET(allGDSl)

Ml <- extractMatrixFromEset(ESETl)

library(RankProd)

# CAUTION: 100 MB, but without huge=TRUE output file is ~ 3gb.
startTime <- proc.time()
res <- RankProducts(data=exprs(ESETl[[3]]),
                    cl = allGDSl[[3]]@dataTable@columns$disease.state %>% as.factor() %>% as.numeric - 1,
                    gene.names = rownames(exprs(ESETl[[3]])),
                    logged = TRUE,
                    na.rm=FALSE,
                    plot=FALSE,
                    rand=123, 
                    huge=TRUE
                    )
endTime <- proc.time()
endTime - startTime


par(pch=20)
par(cex=.1)
par(col="black")
plotRP(res, cutoff=.05)

toped <- topGene(res, logged=TRUE, gene.names = rownames(exprs(ESETl[[3]])), cutoff = .05)


plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '236995_x_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '1557845_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '211122_s_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '242603_x_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '239591_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '1557993_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))

par(mfrow=c(1,1))
plot(exprs(ESETl[[3]])[rownames(exprs(ESETl[[3]])) == '1557993_at'], col=factor(allGDSl[[3]]@dataTable@columns$disease.state))

allGDSl[[3]]@dataTable@table %>% filter( ID_REF == '1557993_at') %>% select(-ID_REF, -IDENTIFIER) %>% c %>% unlist %>% as.numeric %>% log2 %>% plot(., col=factor(allGDSl[[3]]@dataTable@columns$disease.state))

tmp <- as.data.frame(exprs(ESETl[[3]])) %>% rownames_to_column('probe') %>% filter( probe %in% rownames(toped$Table1)[1:30] ) %>% gather(sample, eval, starts_with('GSM')) %>% left_join(., allGDSl[[3]]@dataTable@columns %>% select(sample, disease.state) %>% mutate(sample = as.character(sample) ), by="sample") 
p <- ggplot(tmp) + geom_point(aes(x=sample,y=eval,color=disease.state)) + facet_wrap( ~ probe )
ggsave(filename = paste0(gcl$rootDir, '/caprice/EDA/RankProd-GDS1917-Table1-30.png'), plot=p, width = 16, height=10)

# CAUTION: create ~ 3gb output file.
generes <- RankProducts(data=Ml[[3]],
                    cl = allGDSl[[3]]@dataTable@columns$disease.state %>% as.factor() %>% as.numeric - 1,
                    gene.names = rownames(Ml[[3]]),
                    logged = TRUE,
                    na.rm=FALSE,
                    plot=FALSE,
                    rand=123)

plotRP(generes, cutoff=.05)

toped_gene <- topGene(generes, gene.names = rownames(Ml[[3]]), num.gene = 100)



as.data.frame(Ml[[3]]) %>% add_rownames('UNI') %>% filter( UNI %in% rownames(toped_gene$Table1)[1:10] ) %>% gather(sample, eval, starts_with('GSM')) %>% left_join(., allGDSl[[3]]@dataTable@columns %>% select(sample, disease.state) %>% mutate(sample = as.character(sample) ), by="sample") %>%  ggplot(aes(x=sample,y=eval)) + geom_point(aes(color=disease.state)) + facet_wrap( ~ UNI )

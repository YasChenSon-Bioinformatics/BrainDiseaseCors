library(ranger)

setGlobalConstantList()
loadLibraries()
allGDSl <- download_GDSs(skipv=c('not-GPL570', 'no-disease.state', 'no-control',
                                 'non-binary', 'found-NA'))
ESETl <- convertGDS2ESET(allGDSl)
Ml <- extractMatrixFromEset(ESETl)

univ <- cbind(t(Ml[[1]])) %>% as.data.frame(., stringsAsFactors = FALSE)# %>% mutate(y = allGDSl[[1]]@dataTable@columns$disease.state )
univ$y <- as.factor(allGDSl[[1]]@dataTable@columns$disease.state)

start <- proc.time()
# DO NOT USE FORMULA. protect() error occurred.
rangered3 <- ranger(data = univ[,1:ncol(univ)], num.trees = 10000, mtry = round(sqrt(ncol(univ)),0), min.node.size = 1, scale.permutation.importance = TRUE, importance = "permutation", write.forest = FALSE, save.memory = FALSE, dependent.variable.name = 'y', classification = TRUE, seed=1234)
end <- proc.time()
(end - start)['elapsed']

cands <- names(sort(rangered3$variable.importance, decreasing = TRUE)[1:500])

tst <- cbind(t(Ml[[2]])) %>% as.data.frame(., stringsAsFactors = FALSE)# %>% mutate(y = allGDSl[[1]]@dataTable@columns$disease.state )
tst$y <- as.factor(allGDSl[[2]]@dataTable@columns$disease.state)

rangered4 <- ranger(data = tst[,1:ncol(tst)], num.trees = 10000, mtry = round(sqrt(ncol(tst)),0), min.node.size = 1, scale.permutation.importance = TRUE, importance = "permutation", write.forest = FALSE, save.memory = FALSE, dependent.variable.name = 'y', classification = TRUE, seed=1234)


rangered4$variable.importance

cands4 <- names(sort(rangered4$variable.importance, decreasing = TRUE)[1:500])


# FIXME: invalid test set usage! Need to scale intensities. 

library(randomForest)
x <- randomForest(x=univ[, c('243245_at', '1562245_a_at', '205420_at')], y=factor(allGDSl[[1]]@dataTable@columns$disease.state), importance = TRUE, ntree = 500, )
x <- randomForest(x    =univ[, cands], y    =factor(allGDSl[[1]]@dataTable@columns$disease.state),
                  xtest= tst[, cands], ytest=factor(allGDSl[[2]]@dataTable@columns$disease.state), importance = TRUE, ntree = 20000 )
z <- randomForest(x    = tst[, cands], y    =factor(allGDSl[[2]]@dataTable@columns$disease.state),
                  xtest=univ[, cands], ytest=factor(allGDSl[[1]]@dataTable@columns$disease.state), importance = TRUE, ntree = 20000 )
w <- randomForest(x    = tst[, cands4], y    =factor(allGDSl[[2]]@dataTable@columns$disease.state),
                  xtest=univ[, cands4], ytest=factor(allGDSl[[1]]@dataTable@columns$disease.state), importance = TRUE, ntree = 20000 )

x$inbag


which(names(sort(rangered3$variable.importance, decreasing = TRUE)) == '1563260_at')
which(names(sort(rangered3$variable.importance, decreasing = TRUE)) == '1563260_at')

df <-
left_join(
x$importance %>% as.data.frame %>% rownames_to_column('probe'),
z$importance %>% as.data.frame %>% rownames_to_column('probe'), by="probe") 

df %>% mutate( score = MeanDecreaseGini.x + MeanDecreaseGini.y ) %>% arrange(desc(score))




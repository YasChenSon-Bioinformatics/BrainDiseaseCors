# Reusable functions

download_GDS <- function(GDSnumberv) {
  all_GDS = list()
  j = 1
  for (gds in GDSnumberv){
    message(paste('GDS', gds, sep=''))
    all_GDS[[j]] = getGEO(GEO=paste('GDS', GDSnumberv[j], sep=''), destdir='/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/data')
    j = j + 1
  }
  all_ESET = list()
  j = 1
  for (gds in all_GDS){
    all_ESET[[j]] = assign(paste('ESET', GDSnumberv[j], sep='_'), GDS2eSet(GDS=gds)) 
    j = j + 1
  }
  all_DF = list()
  j = 1;
  for (eset in all_ESET) { 
    m <- exprs(eset)
    pdata <- pData(eset)
    all_DF[[j]] <- assign(paste('df', GDSnumberv[j], sep='_'), cbind(pdata, t(m)))
    j = j + 1
  }
  
  return(all_DF)
}





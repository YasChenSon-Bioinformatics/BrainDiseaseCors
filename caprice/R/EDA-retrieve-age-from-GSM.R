

library(GEOmetadb)
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0046178
GDS4358 <- getGEO('GDS4358', destdir = '/Users/PCUser/Downloads/Rtmp/' )
GDS4358@dataTable@columns

geoConvert( GDS4358@header$dataset_id[1], c('gse','gpl','gsm','gds')[1] )
sfile = getSQLiteFile() # download a copy of the GEOmetadb SQLite file
# con = dbConnect(SQLite(),sfile)
# dbGetQuery(con,'select gds,title,gse,gpl from gds')

queryGEO <- function(
  GSE='GSE35864'
){
  queried <- geoConvert(GSE, 'gsm')
  out <- list()
  for( i in seq_along(queried$gsm$to_acc) ){
    GSM <- getGEO(queried$gsm$to_acc[i], destdir = '/Users/PCUser/Downloads/Rtmp/')
    out[[i]] <- GSM
  }
  a<- sapply(out, function(smpl) {smpl@header$characteristics_ch1})
  colnames(a) <- queried$gsm$to_acc
  GDSmetadf <- t(a)
  return(GDSmetadf)
}

GDS4358metadf <- queryGEO(GDS4358@header$reference_series)


### https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2154

GDS2154 <- getGEO('GDS2154', destdir = '/Users/PCUser/Downloads/Rtmp/' )
GDS2154@dataTable@columns # no age

GDS2154metadf <- queryGEO(GDS2154@header$reference_series)



#GSM876886 <- getGEO("GSM876886", destdir = '/Users/PCUser/Downloads/Rtmp/')
GSM876886@header$characteristics_ch1[2]

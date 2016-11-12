source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()
GDSl <- download_GDSs(skipv=c('not-GPL570', 'blacklist'))

GDS_strv <- sapply(GDSl, function(x){ x@header$dataset_id[1]})

library(GEOmetadb) # not in loadLibraries()

getSQLiteFile()
# If the above function getSQLiteFile() fails, just download from
# https://dl.dropboxusercontent.com/u/51653511/GEOmetadb.sqlite.gz

con <- dbConnect(SQLite(),'/Users/PCUser/Downloads/Rtmp/GEOmetadb.sqlite')
geo_tables <- dbListTables(con)
geo_tables
#dbListFields(conn = con, name = 'gse_gsm')
sapply(geo_tables, function(x) dbListFields(con,x))

dbGetQuery(con, "select * from gds limit 3")
# You can use dbGetQuery here, but you can also access to GEOmetadb.sqlite file by dplyr


# 1 gds <----> 1 gse (GEo Series) <----> n samples <----> n age information

db <- src_sqlite('/Users/PCUser/Downloads/Rtmp/GEOmetadb.sqlite')


GSEv <- tbl(db,'gds') %>% filter( gds %in% GDS_strv ) %>% select(gse) %>% collect() %>% c %>% unlist


partialQuery <- paste0("'", paste(GSEv, collapse="', '"), "'" )
# tbl(db, 'gse_gsm') %>% rename(G=gse) %>% filter( G %in% GSEv )
gsmdf <- dbGetQuery(con, paste("select * from gse_gsm where gse IN(", partialQuery, ')'))

# length(gsmdf$gsm) == length(unique(gsmdf$gsm))

tmp <- tbl(db, 'gsm') %>% filter( gsm %in% gsmdf$gsm) %>% select(gsm, characteristics_ch1) %>% collect


grepl(pattern = 'abc', x = 'jfdsabc')    # grep by Logical
grepl(pattern = 'abc', x = 'jfdsab ')     # grep by Logical

#tmp2 <-
#tmp %>% filter( grepl('[Aa]ge', characteristics_ch1) ) 
agedf <- tmp %>% filter(  (!grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) & grepl('[aA]ge:',characteristics_ch1)) ) %>% rename( age = characteristics_ch1)  %>% mutate( age = str_extract(string=age, pattern = 'age:[ 0-9]*') %>% gsub('[^0-9]', '', .) %>% as.numeric())


hist(agedf$age)
GSE_with_age <- gsmdf %>% filter( gsm %in% agedf$gsm ) %>% .$gse %>% unique
tbl(db,'gds') %>% filter( gse %in% GSE_with_age )
 
SAMPLES_WITH_AGE <- as.character(GDSlist[[4]]@dataTable@columns$sample) # GDS4358
agedf %>% filter( gsm %in% SAMPLES_WITH_AGE )

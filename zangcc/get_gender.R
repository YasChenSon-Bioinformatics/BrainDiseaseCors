source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList(rootDir = '/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree')
loadLibraries()
GDSl <- download_GDSs(skipv=c('not-GPL570', 'blacklist'))

GDS_strv <- sapply(GDSl, function(x){ x@header$dataset_id[1]})

source("https://bioconductor.org/biocLite.R")
biocLite("GEOmetadb")

library(GEOmetadb) # not in loadLibraries()

getSQLiteFile()
# If the above function getSQLiteFile() fails, just download from
# https://dl.dropboxusercontent.com/u/51653511/GEOmetadb.sqlite.gz

con <- dbConnect(SQLite(),'/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree/zangcc/GEOmetadb.sqlite')
geo_tables <- dbListTables(con)
geo_tables
#dbListFields(conn = con, name = 'gse_gsm')
sapply(geo_tables, function(x) dbListFields(con,x))

dbGetQuery(con, "select * from gds limit 3")
# You can use dbGetQuery here, but you can also access to GEOmetadb.sqlite file by dplyr


# 1 gds <----> 1 gse (GEo Series) <----> n samples <----> n age information

db <- src_sqlite('/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree/zangcc/GEOmetadb.sqlite')


GSEv <- tbl(db,'gds') %>% filter( gds %in% GDS_strv ) %>% select(gse) %>% collect() %>% c %>% unlist


#Do some formating like adding quote for characters
partialQuery <- paste0("'", paste(GSEv, collapse="', '"), "'" )
# tbl(db, 'gse_gsm') %>% rename(G=gse) %>% filter( G %in% GSEv )
gsmdf <- dbGetQuery(con, paste("select * from gse_gsm where gse IN(", partialQuery, ')'))

# length(gsmdf$gsm) == length(unique(gsmdf$gsm))

tmp <- tbl(db, 'gsm') %>% filter( gsm %in% gsmdf$gsm) %>% select(gsm, characteristics_ch1) %>% collect


grepl(pattern = 'abc', x = 'jfdsabc')    # grep by Logical
grepl(pattern = 'abc', x = 'jfdsab ')     # grep by Logical

genderdf <-
  tmp %>% filter( 
    (!grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) &
       grepl('([sS]ex|[gG]ender)', characteristics_ch1))) %>%
  rename( gender = characteristics_ch1) %>%
  mutate( gender = str_extract(string=gender, pattern = '([sS]ex|[gG]ender): ([mM]|[mM]ale|[fF]|[fF]emale)') ) %>%
  mutate( gender = gsub('.* ','', gender) %>% toupper )

# table(genderdf$gender)

GSE_with_gender <- gsmdf %>% filter( gsm %in% genderdf$gsm ) %>% .$gse %>% unique
tbl(db,'gds') %>% filter( gse %in% GSE_with_gender )


SAMPLES_WITH_GENDER_1 <- as.character(GDSl[[6]]@dataTable@columns$sample) #GDS4136 
SAMPLES_WITH_GENDER_2 <- as.character(GDSl[[4]]@dataTable@columns$sample) #GDS4358
GDS4136_join_table <- genderdf %>% filter( gsm %in% SAMPLES_WITH_GENDER_1 )
GDS4358_join_table <- genderdf %>% filter( gsm %in% SAMPLES_WITH_GENDER_2 )

after_join_GDS4136 <- left_join(GDS4136_join_table, GDSl[[6]]@dataTable@columns, by = c('gsm' = 'sample')) %>% mutate( gender = gsub('.* ','', gender) %>% toupper) #left join GDS4136

after_join_GDS4358 <- left_join(GDS4358_join_table, GDSl[[4]]@dataTable@columns, by = c('gsm' = 'sample')) %>% mutate( gender = gsub('.* ','', gender) %>% toupper) #left join GDS4136 


# Consideration of gender
ESET4136 <- GDS2eSet(GDSl[[6]], do.log2 = TRUE)
MATRIX4136 <- as.matrix( ESET4136 )
gene_expression_values <- exprs(ESET4136)
dz <- GDSl[[6]]@dataTable@columns$disease.state
gender <- after_join_GDS4136$gender
lmfitted2 <- lmFit(MATRIX4136, design = model.matrix( ~  dz + gender))
ebayesd2 <- eBayes(lmfitted2)                
topped2 <- topTable(ebayesd2, number = 400)

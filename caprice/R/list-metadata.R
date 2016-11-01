
function(){
  
  source('~/BrainDiseaseCors/caprice/R/analyze-GPL570.R')
  setGlobalConstantList()
  loadLibraries()
  GDSl <- download_GDSs(skipv=c('not-GPL570'))
  ESETl <- convertGDS2ESET(GDSl)
  Ml <- extractMatrixFromEset(ESETl)
  
}

glimpse(this@dataTable@columns)

queryGEO <- function(
  header
){
  queried <- geoConvert(header$reference_series, 'gsm')
  out <- list()
  entry <- 
  for( i in seq_along(queried$gsm$to_acc) ){
    GSM <- getGEO(queried$gsm$to_acc[i], destdir = '/Users/PCUser/Downloads/Rtmp/')
    out[[i]] <- GSM
  }
  a<- sapply(out, function(smpl) {smpl@header$characteristics_ch1})
  
  if ( length(dim(a)) < 2 ){
    names(a) <- queried$gsm$to_acc
    GDSmetadf <- a
  } else {
    colnames(a) <- queried$gsm$to_acc
    GDSmetadf <- t(a)    
  }

  return(list(gds = header$dataset_id[1], title = header$title,
              desc = header$description, df =GDSmetadf, item = out[[1]]@header$characteristics_ch1))
}

metal <- list(); for ( i in seq_along(GDSl) ){
  metal[[i]] <- queryGEO(GDSl[[i]]@header)
}

message(paste0(sapply(metal, function(x) paste0(x$gds, '\t', substr(x$title,1,35), '\t', paste0(x$item[grepl('age|gender|[Ss]ex|male',x$item)], collapse='\t') ) ),collapse='\n'))
#    age : 9 out of 17
# gender:  9 out of 17

library(GEOmetadb)

con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
geo_tables <- dbListTables(con)
geo_tables
sapply(geo_tables, function(x) dbListFields(con,x))

trials_for_inspection <- function(){
  dbGetQuery(con,'PRAGMA TABLE_INFO(gpl)')
  rs <- dbGetQuery(con,'select gsm, gpl, characteristics_ch1 from gsm where gpl = "GPL570" and characteristics_ch1  like "%age=%"  limit 100')
  dbGetQuery(con,'select count(*) from gsm where gpl = "GPL570" and (characteristics_ch1  like "%age=%" or characteristics_ch1 like "age :")')
  rs
  
  
  
  dbGetQuery(con,'select count(*), sum(sample_count) from gds where gpl = "GPL570"')                                # 527 GDSs, 13381 samples
  dbGetQuery(con,'select count(*), sum(sample_count) from gds where gpl = "GPL570" and value_type != "log2 ratio"') # 526 GDSs, 13347 samples
  dbGetQuery(con,'select distinct(value_type) from gds where gpl = "GPL570"') # 527 GDSs, 13381 samples
  
  # GDS1411 refers to GSE2555. GSE2555 contains TWO platforms (GPL97 and GPL570). Thus we need to filter out some samples in GSE2555
  GSM_in_GDS_on_GPL570 <- # 
    dbGetQuery(con,
               paste('SELECT B.gsm, B.gse, B.gds, B.sample_count  ',
                     'FROM ((SELECT AA.gsm, AA.gse, A.gds, A.sample_count FROM ((SELECT gse, gds, sample_count FROM gds WHERE gpl = "GPL570") AS A JOIN gse_gsm USING (gse)) AS AA)',
                     '      JOIN gsm USING (gsm)',
                     '      ) AS B WHERE gpl = "GPL570"'))        
  dbGetQuery(con, 'SELECT gds, count(*) FROM gds GROUP BY gds ')
  samples_in_gds <- dbGetQuery(con,'select , sum(sample_count) from gds where gpl = "GPL570"')                                # 526 GDSs, 13347 samples
  
  library(dplyr)
  db = src_sqlite('GEOmetadb.sqlite')
  tbl(db,'gsm') %>% filter( gpl == 'GPL570') %>% select(gsm, characteristics_ch1) %>% filter( characteristics_ch1 %like% "%age=%" )
  tbl(db,'gsm') %>% filter( gpl == 'GPL570') %>% select(gsm, characteristics_ch1) %>% filter(!regexp("^1", characteristics_ch1))
  tbl(db,'gsm') %>% filter( gpl == 'GPL570') %>% select(gsm, characteristics_ch1) %>% filter( characteristics_ch1 %like% "%age:%" | characteristics_ch1 %like% "%Age:%" ) %>% head(1000) %>% View
  
  tbl(db,'gsm') %>% filter( gpl == 'GPL570') %>% select(gsm, characteristics_ch1) %>% filter( characteristics_ch1 %like% "%age:%" | characteristics_ch1 %like% "%Age:%" ) %>% collect %>% nrow # 27589
  tbl(db,'gsm') %>% filter( gpl == 'GPL570') %>% select(gsm, characteristics_ch1) %>% filter( characteristics_ch1 %like% "%gender:%" | characteristics_ch1 %like% "%Gender:%" | characteristics_ch1 %like% "%male%" ) %>% collect %>% nrow # 25963
  
}

library(pbapply)

GDS_on_GPL570        <- dbGetQuery(con,'select gds, gse, sample_count from gds where gpl = "GPL570"')
GSM_on_GPL570_in_GDS <- pbapply::pbsapply(GDS_on_GPL570$gse, FUN = function(gse) geoConvert(gse, 'gsm')) # takes >10 minutes
save(GSM_on_GPL570_in_GDS, file='BrainDiseaseCors/caprice/DATA/GSM_on_GPL570_in_GDS.Rdata') # 1.5 MB

GSM_duplicated <- sapply(GSM_on_GPL570_in_GDS, function(df) df$to_acc ) %>% unlist # 15222 samples
GSM_unique     <- unique(GSM_duplicated)                                           # 13979 samples ... GDS seems remove some of them. But ignore them at this moment

tbl(db,'gsm') %>% filter( gsm %in% GSM_unique ) %>% select(gsm, characteristics_ch1) %>% filter( characteristics_ch1 %like% "%age:%" | characteristics_ch1 %like% "%Age:%" )
# GDS3341 has Tumor st'age'. So ' age' might be better
# GDS5029: ine'age': colorectal;	cell line: SW480;	treatment: none
# GDS3936
gds2gsm <- 
left_join(
  tbl(db,'gsm') %>% filter( gsm %in% GSM_unique ) %>% select(gsm, characteristics_ch1) %>%
    filter( characteristics_ch1 %like% "% age:%" |
            characteristics_ch1 %like% "%Age:%" |  # GDS2822 has 'maternalAge:'
            characteristics_ch1 %like% "%\tAge:%" |
            characteristics_ch1 %like% "%\tAge:%" 
    ) %>% collect,
  do.call(rbind, GSM_on_GPL570_in_GDS) %>% rename( gsm = to_acc, gse = from_acc ),
  by='gsm'
) %>% left_join(., GDS_on_GPL570, by='gse')
gds2gsm %>% group_by(gds) %>% filter( !(grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) & !grepl(' age|Age:',characteristics_ch1)) )  %>% slice(1) # 63 GDSs have 'Age' columns. ~50 has gender columns
GDS_with_age_metadata <- gds2gsm %>% group_by(gds) %>% filter( !(grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) & !grepl(' age|Age:',characteristics_ch1)) )  %>% slice(1)

ageGDSl<- download_GDSs( GDSnumberv = gsub('GDS', '', GDS_with_age_metadata$gds), skipv = c('not-GPL570'), suppress = FALSE ) # 4GB!
sapply(ageGDSl, function(x) 'age' %in% colnames(x@dataTable@columns) ) %>% sum()
sum(sapply(ageGDSl, function(x) ifelse('age' %in% colnames(x@dataTable@columns),0,as.numeric(x@header$sample_count)) ))
sapply(ageGDSl, function(x) ifelse('age' %in% colnames(x@dataTable@columns), '', x@header$dataset_id[1]))
tmp <- sapply(ageGDSl, function(x) ifelse('age' %in% colnames(x@dataTable@columns), '', x@header$dataset_id[1])) %>% gsub("GDS","",.) 
paste(tmp[nchar(tmp)!=0],collapse=', ')

GDS_on_GPL570_with_no_age_in_original_experimental_settings_but_has_info_in_corresponding_GSE <- c(1344, 1665, 1989, 2052, 2154, 2416, 2418, 2453, 2486, 2491, 2495, 2534, 2609, 2628, 2652, 2822, 2887, 2970, 3256, 3423, 3463, 3496, 3554, 3630, 3688, 3841, 3880, 3901, 3954, 3955, 3961, 4053, 4102, 4147, 4195, 4279, 4282, 4296, 4297, 4327, 4345, 4354, 4358, 4467, 4468, 4473, 4477, 4519, 4532, 4580, 4620, 4837, 4854, 4901) # 54

gds2gsm %>% group_by(gds) %>% filter( !(grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) & !grepl(' age|Age:',characteristics_ch1)) )  %>% slice(1) %>% filter(gds %in% sapply(ageGDSl, function(x) ifelse('age' %in% colnames(x@dataTable@columns), '', x@header$dataset_id[1])))

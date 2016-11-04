
GDS_on_GPL570_with_no_age_in_original_experimental_settings_but_has_info_in_corresponding_GSE <-
  c(1344, 1665, 1989, 2052, 2154, 2416, 2418, 2453, 2486, 2491,
    2495, 2534, 2609, 2628, 2652, 2822, 2887, 2970, 3256, 3423,
    3463, 3496, 3554, 3630, 3688, 3841, 3880, 3901, 3954, 3955,
    3961, 4053, 4102, 4147, 4195, 4279, 4282, 4296, 4297, 4327,
    4345, 4354, 4358, 4467, 4468, 4473, 4477, 4519, 4532, 4580,
    4620, 4837, 4854, 4901) # 54
GDS_with_ages <- paste0('GDS', GDS_on_GPL570_with_no_age_in_original_experimental_settings_but_has_info_in_corresponding_GSE)
# GEOmetadb is downloaded from https://dl.dropboxusercontent.com/u/51653511/GEOmetadb.sqlite.gz
library(GEOmetadb)
library(pbapply)
dbpath <- '~/BrainDiseaseCors/caprice/DATA/GEOmetadb.sqlite'
con <- dbConnect(SQLite(),'~/BrainDiseaseCors/GEOmetadb.sqlite') # assume GEOmetadb.sqlite is extracted and placed in home folder (~)
GSE_ages <- pbapply::pbsapply(GDS_with_ages, function(gds) geoConvert(gds, 'gse', sqlite_db_name = dbpath))
GSEages <- sapply(GSE_ages, function(x) x$to_acc)
GSM_on_GPL570_in_GDS <- pbapply::pbsapply(GSEages, FUN = function(gse) geoConvert(gse, 'gsm', sqlite_db_name = dbpath)) # takes >3 minutes
#save(GSM_on_GPL570_in_GDS, file='BrainDiseaseCors/caprice/DATA/GSM_with_ages.Rdata') # 1.5 MB

GSM_duplicated <- sapply(GSM_on_GPL570_in_GDS, function(df) df$to_acc ) %>% unlist # 1587 samples
GSM_unique     <- unique(GSM_duplicated)                                           # 1514 samples

db = src_sqlite('~/BrainDiseaseCors/caprice/DATA/GEOmetadb.sqlite')

gsmdf <-
  tbl(db,'gsm') %>% filter( gsm %in% GSM_unique ) %>% dplyr::select(gsm, characteristics_ch1) %>%
  filter( characteristics_ch1 %like% "% age:%"  |
            characteristics_ch1 %like% "% Age:%"  |
            characteristics_ch1 %like% "%age:%"   |  # GDS2822 has 'maternalAge:'
            characteristics_ch1 %like% "%Age:%"   |  # GDS2822 has 'maternalAge:'
            characteristics_ch1 %like% "%\tAge:%" |
            characteristics_ch1 %like% "%\tAge:%" 
  ) %>% collect # 1418 GSM

cands <- gsmdf %>% dplyr::rename(m=characteristics_ch1) %>% dplyr::filter( !(grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',m) & !grepl(' age|Age:',m)) ) %>%
  mutate(age = str_extract(m,'[Aa]ge[: 0-9]+') %>% gsub('[^0-9]','',.)) %>%
  filter(age != '') %>% # exclude non-numeric records such as 'age: Adult' (GSM656466) or 'Age: < 35 yrs.' (GSM304262)
  mutate(age = ifelse(age==358, 35, as.numeric(age) )) # GSM277455 has incorrect age records


# cands_GSE <- pbapply::pbsapply(cands$gsm, FUN = function(gsm) geoConvert(gsm, 'gse', sqlite_db_name = dbpath)) # too slow
# geoConvert('GSM135245', 'gse', sqlite_db_name = dbpath)
left_join(do.call(rbind,GSM_on_GPL570_in_GDS) %>% dplyr::rename( gse = from_acc, gsm = to_acc), 
          do.call(rbind,GSE_ages) %>% group_by(to_acc) %>% dplyr::rename(gse = to_acc, gds = from_acc),  # summarize( GDS = paste0(from_acc, collapse=',')) %>% 
          by='gse') %>%
  left_join(cands, ., by='gsm') %>% .$gds %>% unique

write.csv(cands[, c('gsm','age')], '~/BrainDiseaseCors/caprice/MAP/gsm2age.csv', quote=FALSE, row.names=FALSE)

GDS_with_age_metadata <- 
  c(1344, 1665, 1989, 2052, 2154, 2416, 2418, 2453, 2486, 2491,
    2495, 2534, 2609, 2628, 2652, 2822, 2887, 2970, 3423, 3463,
    3496, 3554, 3630, 3688, 3841, 3880, 3954, 3955, 3961, 4053,
    4102, 4147, 4195, 4279, 4282, 4296, 4327, 4345, 4354, 4358,
    4467, 4468, 4473, 4477, 4519, 4532, 4580, 4620, 4837, 4854, 4901)
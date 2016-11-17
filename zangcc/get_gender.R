source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()
GDSl <- download_GDSs(skipv=c('not-GPL570', 'blacklist'))

GDS_strv <- sapply(GDSl, function(x){ x@header$dataset_id[1]})

source("https://bioconductor.org/biocLite.R")
biocLite("GEOmetadb")

library(GEOmetadb) # not in loadLibraries()

getSQLiteFile()
# If the above function getSQLiteFile() fails, just download from
# https://dl.dropboxusercontent.com/u/51653511/GEOmetadb.sqlite.gz

if ( Sys.info()['user'] == 'PCUser' ) {
    metadb_path <- '/Users/PCUser/Downloads/Rtmp/GEOmetadb.sqlite'
} else {
    metadb_path <- '/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree/zangcc/GEOmetadb.sqlite'
}
con <- dbConnect(SQLite(), metadb_path)
geo_tables <- dbListTables(con)
geo_tables
#dbListFields(conn = con, name = 'gse_gsm')
sapply(geo_tables, function(x) dbListFields(con,x))

dbGetQuery(con, "select * from gds limit 3")
# You can use dbGetQuery here, but you can also access to GEOmetadb.sqlite file by dplyr


# 1 gds <----> 1 gse (GEo Series) <----> n samples <----> n age information

db <- src_sqlite(metadb_path)


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

i_GDS4136 <- which(sapply(GDSl, function(thisGDS) thisGDS@header$dataset_id[1] ) == 'GDS4136' ) # 6
i_GDS4358 <- which(sapply(GDSl, function(thisGDS) thisGDS@header$dataset_id[1] ) == 'GDS4358' ) # 4

GDS4136columns <- GDSl[[ i_GDS4136 ]]@dataTable@columns %>% mutate( sample = as.character(sample) )
GDS4358columns <- GDSl[[ i_GDS4358 ]]@dataTable@columns %>% mutate( sample = as.character(sample) )

GDS4136_join_table <- genderdf %>% filter( gsm %in% GDS4136columns$sample )
GDS4358_join_table <- genderdf %>% filter( gsm %in% GDS4358columns$sample )

after_join_GDS4136 <- left_join(GDS4136_join_table, GDS4136columns, by = c('gsm' = 'sample')) #left join GDS4136
after_join_GDS4358 <- left_join(GDS4358_join_table, GDS4358columns, by = c('gsm' = 'sample')) #left join GDS4136 


# Consideration of gender
ESET4136 <- GDS2eSet(GDSl[[ i_GDS4136 ]], do.log2 = TRUE)
MATRIX4136 <- as.matrix( ESET4136 )
gene_expression_values <- exprs(ESET4136)
dz <- GDS4136columns$disease.state
gender <- after_join_GDS4136$gender
lmfitted2 <- lmFit(MATRIX4136, design = model.matrix( ~  dz + gender ))
ebayesd2 <- eBayes(lmfitted2)                
topped2 <- topTable(ebayesd2, number = 400)

# take 1
top_probev <- rownames(topped2)[1:20]
# take 2
top_probev <- names(sort(abs(lmfitted2$coefficients[, 'M']), decreasing = TRUE)[21:40])
# take 3
malev <- after_join_GDS4136$gender == 'M'
tmp <- rowMeans(gene_expression_values[, malev]) - rowMeans(gene_expression_values[, ! malev])
top_probev <- names(sort(abs(tmp), decreasing = TRUE)[1:20])

# diagnostics
png('zangcc/img/GDS4136-top20probe.png', width = 960, height = 960 ) # width and height are pixel
    gene_expression_values[top_probev, ] %>%
    as.data.frame %>% # dplyr cannot handle a class 'matrix'
    rownames_to_column('probe') %>%
    gather(key=smpl, value=eval, -probe) %>% # See: https://blog.rstudio.org/2014/07/22/introducing-tidyr/
    left_join(. , after_join_GDS4136, c('smpl' = 'gsm')) %>% # we need disease.state later
    ggplot() +
    geom_text(aes(x=smpl, y=eval, label=gender, color=disease.state)) +
    facet_wrap( ~ probe, nrow =5) + # create subplots for each probe
    theme( axis.text.x = element_text(angle=90) ) # rotate sample label in x axis
dev.off() # finish outputting to the above png file

colnames(lmfitted2$coefficients) <- gsub('\\(|\\)|^dz|gender','',colnames(lmfitted2$coefficients))

hypothesis <- lmfitted2$coefficients[top_probev, ] %>% as.data.frame %>%
    rownames_to_column('probe') %>% rename( control = Intercept ) %>%
    mutate_each(funs( . = . + control ), ends_with('_stage')) %>%
    gather(key=dz_state, value=predicted_eval, -M, -probe)

sample_indexes <-
    table(GDS4136columns$disease.state) %>%
    cumsum %>% as.data.frame %>%   # cumulative sum for later visualization
    rename_( xend = '.' ) %>%  # dplyr::rename_() can handle 'quoted' variable names, while dplyr::rename() cannot
    rownames_to_column('dz_state') %>%
    mutate( x = lag(xend, n=1, default = 0) + 1)

vised_hypos <- left_join(hypothesis, sample_indexes, by='dz_state')  # Visualized Hypotheses


png('zangcc/img/GDS4136-top20probe-gender.png', width = 1440, height = 1440 ) # width and height are pixel
gene_expression_values[ top_probev, ] %>%
    as.data.frame %>% # dplyr cannot handle a class 'matrix'
    rownames_to_column('probe') %>%
    gather(key=smpl, value=eval, -probe) %>% # See: https://blog.rstudio.org/2014/07/22/introducing-tidyr/
    left_join(. , after_join_GDS4136, c('smpl' = 'gsm')) %>% # we need disease.state later
    ggplot() +
    geom_text(aes(x=smpl, y=eval, label=gender, color=disease.state)) +
    geom_segment(aes(x= x        ,xend=   xend   , y=predicted_eval,yend=predicted_eval    , color=dz_state), data=vised_hypos) + # For Male
    geom_segment(aes(x=(x+xend)/2,xend=(x+xend)/2, y=predicted_eval,yend=predicted_eval + M, color=dz_state), data=vised_hypos,
                 arrow = arrow(length = unit(0.2,"cm") ) ) + # For Female
    facet_wrap( ~ probe, nrow =5) + # create subplots for each probe
    theme( axis.text.x = element_text(angle=90) ) + # rotate sample label in x axis
    labs(title='GDS4136 topTable() highest 20 probe expression values\nHorizontal Bars denote Hypotheses. Arrows for Female effect')
dev.off() # finish outputting to the above png file


names(sort(abs(lmfitted2$coefficients[, 'M']), decreasing = TRUE)[1:20])



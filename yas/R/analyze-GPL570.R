all_GDSv <- c(5204,4879,4859,4838,4758,4532,4522,4523,4477,4414,4358,4231,4218,4154,4136,4135,3834,3502,3459,3345,3129,3128,3113,3110,3069,2978,2941,2821,2795,2613,2191,2190,2154,1962,1917,1912,1835,1816,1815,1813,1726,1253,1096,1085,910,909,833,810,707,596,564,426,232,181)

pfc = c(2190, 3502, 4414, 4523) #, 4532) - not in same region
options(max.print=1000)
gcl <- list(
    affy2uni = '/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/yas/MAP/affy2uni.txt'
)

library(GEOquery)
library(limma)   # as.matrix.ExpressionSet is defined in limma

# Download all the datasets
all_GDS <- list()
i <- 1
for ( no in all_GDSv ) {
    if (Sys.info()['user'] == 'PCUser'){
        # To avoid downloading the same data again and again.
        # Sometimes downloaded data were corrupted without any warning. Scary
        thisGDS <- getGEO(GEO=paste0("GDS", no), destdir = '/Users/PCUser/Downloads/Rtmp')
    } else {
        thisGDS <- getGEO(GEO=paste0("GDS", no))
    }
    #Sys.sleep(5)
    if( ! thisGDS@header$platform %in% c('GPL570') ) {
        message('Skip GDS ', no, ' because of platform')
        next
    }
    if( ! 'disease.state' %in% colnames(thisGDS@dataTable@columns) ){
        message('Skip GDS ', no, ' because of no disease.state')
        next
    }
    all_GDS[[i]] <- thisGDS
    i <- i + 1
}

# Note: some rows like 1320_at or 1405_i_at in GDS4522 are <NA>.
# 

GDS_on_GPL570v <- sapply( all_GDS, function(gds) {gds@header$dataset_id[1]} )

# all GDS in all_GDSv               ~ 2.5 GB
# only GDS with the platform GPL570 ~ 1.2 GB

# Extract all to ESETs
all_ESET <- list()
for (i in seq_along(all_GDS)) {
    all_ESET[[i]] <- GDS2eSet(GDS=all_GDS[[i]], do.log2 = FALSE)
    
    # Note: log2(x) is undefined and return NaN if x is negative.
}
# all_ESET in all_GDSv              ~ 1.6 GB
# only GDS with the platform GPL570 ~ 0.7 GB 

renamed_colname <- 'affy'
affy2uni_df <-
    read.delim(gcl$affy2uni, sep='\t') %>%          # %>% is a pipe. Output of the left is passed to the right
    dplyr::select(-X) %>%                           # avoid confliction between MASS::select() and dplyr::select()
    separate_rows(UniProt.Accession, sep = ';') %>% # separate multiple values like 'Q5TK75; Q6IUU8; Q6V4Z6;'.
    filter( UniProt.Accession != '-' ) %>%          # NB: some Affy ID has no corresponding Uniplot IDs. Strange.
    # This means that our following results are just within some portion of all genes.
    rename( uni = UniProt.Accession ) %>%
    rename_(.dots = setNames('Affy.ID', renamed_colname)) # by using rename_() not rename(), we can use string variables


newAttrName <- 'GDS'

all_M <- list()
i <- 5 # something bad is happening in this dataset, GDS4231
for (i in seq_along(all_ESET)) { 
    tmp <- as.matrix(all_ESET[[i]]) # NB: call library(limma); and library(affy) before
    merged <- merge(tmp, affy2uni_df, by.x = 'row.names', by.y = renamed_colname, all=FALSE)
    colnames(merged)[colnames(merged) == 'Row.names'] <- renamed_colname  
    mapped <-
        merged %>%
        dplyr::select(-starts_with(renamed_colname)) %>% # remove unnecessary affymetrix ID
        dplyr::group_by(uni) %>%                         # the next operation is done for each Uniplot ID group
        dplyr::summarize_each(funs(mean))                # calculate arithmetic mean for each column
    
    all_M[[i]] <- as.matrix( mapped %>% dplyr::select(-uni) )  # e.g. [ 65107 genes x 41 SaMples ]
    rownames(all_M[[i]]) <- mapped$uni
    attr(all_M[[i]], newAttrName) <- all_ESET[[i]]@experimentData@other$dataset_id[1]
    
    #all_DF[[i]] <- as.data.frame(all_ESET[[i]]) 
    #
    # Applying as.data.frame is not a good idea.
    # It's because R doesn't allow to use rownames starting with numbers.
    # '1007_s_at' or '121_at' is automatically converted to 'X1007_s_at' or 'X121_at'.
    #
    # It induces some errors
    # '[E]specially when you've been looking at your screen for 20 straight hours' :)
    #
    # You can check by
    # m <- as.matrix(    all_ESET[[1]])  # NB: call library(limma); and library(affy) before
    # d <- as.data.frame(all_ESET[[1]])
    # data.frame( matrix_name = rownames(m)[1:5], df_name = colnames(d)[1:5] )
    # 
    # I think this is the reason TA uses as.matrix in place of as.data.frame
}
# all_M in all_GDSv               ~ 600 MB
# all_M with the platform GPL570  ~ 300 MB

sapply(all_M, dim)
sapply(all_M, function(x) {sum(is.na(x))} )

######################################################################
# Notes:

# Certain GDS's have        Funny column names
# 1085, 1813, 3834      "x1", "x2", etc --> what probes are these?
# 1835                  x882, x880
# 3113                  x156427, x139282
# 4758, 4879            X7896736
# 833                   175019765, 175019766
# 909, 910              x559, x558 


# Limit our analyses to GPL570 datasets (most common, and most comprehensive array)
#gpl570 = list(df_5204,df_4838,df_4532,df_4523,df_4522,df_4477,df_4358,df_4231,df_4218,df_4136,df_4135,df_3502,df_2821,df_2795,df_2154,df_1962,df_1917)
# Done in above


# Find NAs, NANs (NANs are a todo)
nan_count <-function (x) sapply(x, function(y) sum(is.nan(y)))
i <- 1 
for (i in seq_along(all_M) ) {
    x <- all_M[[i]]
    message(sprintf("%s %6d x %3d [ %5s ] ", attr(x, newAttrName), nrow(x), ncol(x), any(is.na(x))))
    # print(any(is.nan(df)))
}
# why some entries are lacked? like

#####################################################################
# Notes:

# gpl570 dataset        variable count        NAs?        NANs?
# df_5204                   54679               n           
# df_4838                   54678               n
# df_4532                   54680               n
# df_4523                   45680               n
# df_4522                   54680               y
# df_4477                   54678               y
# df_4358                   54679               n
# df_4231                   54679               y
# df_4218                   54680               y
# df_4136                   54679               n
# df_4135                   54679               y
# df_3502                   54681               y
# df_2821                   54679               y
# df_2795                   54679               n
# df_2154                   54678               n
# df_1962                   54679               y
# df_1917                   54678               n

# Find the probes common to all df's

# common_probe_ids = rownames(all_M[[1]])
# for (m in all_M) {
#     common_probe_ids <- intersect(common_probe_ids, rownames(m))
#     message(length(common_probe_ids))
# }
# affyIDs <- lapply( all_GDS, function(x) { x@dataTable@table$ID_REF} ) %>% unlist %>% unique
# clip <- pipe("pbcopy", "w")                       
# write.table(paste(affyIDs, collapse='\n'), file=clip)                               
# close(clip)

# NB: 54677 Common Columns when applied as.data.frame,
#     which contains 'sample' and 'description' as last 2 columns
# When applied to matrix, 54675 probe ids are common. 


# Subset datasets to just disease states

library(dplyr)
library(tidyr)
# experimental condition vector
exp_condv <- sapply( all_GDS, function(gds) { colnames(gds@dataTable@columns) }) %>% unlist %>% unique 
cond_df <- matrix(data = 0, # all-zero amtrix
                  nrow = length(all_GDS),
                  ncol = length(exp_condv),
                  dimnames = list(GDS_on_GPL570v, exp_condv)) %>% data.frame
for ( i in seq_along(all_GDS) ){
    exp_cond_of_this_dataset <- colnames(all_GDS[[i]]@dataTable@columns)
    cond_df[i, ] <- as.numeric(exp_condv %in% exp_cond_of_this_dataset)
}
cdf <- cond_df %>% add_rownames(var='GDS') %>% gather(key = meta, value = bit, -starts_with('GDS'))
experimental_bitmap <-
    left_join(cond_df %>% add_rownames(var='GDS') %>% gather(key = meta, value = bit, -matches('GDS')),
              cdf %>% group_by(meta) %>% summarize( nentry = sum(bit) ), by="meta") %>%
    mutate(meta=paste(nentry,meta)) %>% ggplot(aes(x=gsub('GDS','',GDS), y=reorder(meta,nentry))) +
    geom_tile(aes(fill=as.factor(bit))) + scale_fill_grey()
print(experimental_bitmap)

selectConditionv <- function( m, g = 'GDS5204' ) {
    condition <- list()
    if ( g %in% paste0('GDS', c(1917, 1962, 2154, 2795, 2821, 3502, 4135, 4136, 4218, 4231,
                                4358, 4522, 4523, 4838)) ) {
        condition$col <- 'disease.state'
        condition$val <- 'control'
    } else if (g == 'GDS5204') { # aging / gender
        condition <- NULL
    } else if (g == 'GDS4477') { # genotype variation
        condition <- NULL
    } else if (g == 'GDS4532') { # tissue / development stage
        condition <- NULL
    }
    return(condition)
}

# Compute correlations within samples from each dataset
disease_ms = list()
i <- k <- 2
for (k in seq_along(all_M)) {
    
    m <- all_M[[k]]
    g <- all_GDS[[k]]
    # Subset to just DZ instances
    condition <- selectConditionv( m, attr(m, 'GDS') )
    
    message(attr(m, 'GDS'), appendLF = FALSE)
    if ( is.null(condition) ){
        message('  SKIP ')
        next  
    }else {
        message()
    }
}

disease_samples <-
    g@dataTable@columns[ g@dataTable@columns[, condition$col ] != condition$val, ]$sample

# NB: some GDSs do not have disease.state like GDS5204 (aging)
disease_ms[[i]] = m[ , colnames(m) %in% disease_samples]
i <- i + 1
}



# Compute correlations within samples from each dataset
cors = list()
for (i in seq_along(disease_ms) ) {
    message(i, "-----------------------------------------------")
    # Compute corr
    corre = cor(disease_ms[[i]])
    cors[[i]] = corre
}

names(cors) = c('df_5204','df_4838','df_4532','df_4523','df_4522','df_4477','df_4358','df_4231','df_4218','df_4136','df_4135','df_3502','df_2821','df_2795','df_2154','df_1962','df_1917')
attach(cors)
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
plot(cors[[2]], main="asdf")
plot(cors[[4]], main="asdf")
plot(cors[[7]], main="asdf")
plot(cors[[10]], main="asdf")
plot(cors[[14]], main="asdf")
plot(cors[[15]], main="asdf")
plot(cors[[17]], main="asdf")



# Average across rows for each disease state
avgs = list()
i = 1
for (df in disease_dfs) {
    print(i)
    print(dim(df))
    print("-----------------------------------------------")
    colMeans = colMeans(df, na.rm=TRUE)
    avgs[[i]] = colMeans
    i = i + 1
}



require(reshape2)
worked = list(4,5,7,8,9,10,11,12,13,14,15,16,17) # (also 2 worked)
worked_2 = list(9,11,14,15,16)
avg_dz_exp_vals = as.data.frame(avgs[2])
for (i in worked_2) {
    avg_dz_exp_vals <- cbind(avg_dz_exp_vals, as.data.frame(avgs[i]))
}



install.packages("corrplot")
library(corrplot)
# Compute correlation between diseases expression values
c = cor(avg_dz_exp_vals, use="complete.obs")
names = c("GDS4838 - Neuroectodermal Tumors",
          "GDS4218 - Multiple Sclerosis Brain Lesions",
          "GDS4135 - Temporal Cortex, Ageing Study",
          "GDS2795 - Alzheimer's Disease",
          "GDS2154 - Inflammatory Dilated Cardiomyopathy",
          "GDS1962 - Glioma tissue")
colnames(c) <- rownames(c) <- names
corrplot.mixed(c, t1.pos="r", t1.col="blue", c1.srt=60, c1.pos="r", cl.align.text="r", mar=c(1,1,1,1), height=1600, width=1600)




pData(ESET_3502)
exprs(ESET_3502)
assayDataElement(ESET_3502)


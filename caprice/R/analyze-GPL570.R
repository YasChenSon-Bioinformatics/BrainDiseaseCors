all_GDSv <- c(5204,4879,4859,4838,4758,4532,4522,4523,4477,4414,4358,4231,4218,4154,4136,4135,3834,3502,3459,3345,3129,3128,3113,3110,3069,2978,2941,2821,2795,2613,2191,2190,2154,1962,1917,1912,1835,1816,1815,1813,1726,1253,1096,1085,910,909,833,810,707,596,564,426,232,181)

pfc = c(2190, 3502, 4414, 4523) #, 4532) - not in same region
options(max.print=1000)
gcl <- list(
    affy2uni = '/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/MAP/affy2uni.txt'
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
        message('     Skip GDS ', no, ' because of platform')
        next
    }
    if( ! 'disease.state' %in% colnames(thisGDS@dataTable@columns) ){
        message('     Skip GDS ', no, ' because of no disease.state')
        next
    }
    thisGDS@dataTable@columns$disease.state <- gsub(' ','_',thisGDS@dataTable@columns$disease.state)
    if( ! any(grepl('control', thisGDS@dataTable@columns$disease.state)) ){
        message('     Skip GDS ', no, ' because of no control')
        #next
    }
    if( length(unique(thisGDS@dataTable@columns$disease.state)) != 2 ){
        message('     Skip GDS ', no, ' because of control values')
        #next
    }
    if( any(is.na(thisGDS@dataTable@table)) ) {
        message('     Skip GDS ', no, ' because of NA')
        #next
    }
    all_GDS[[i]] <- thisGDS
    i <- i + 1
}

function(){
    tmpdf <- sapply(all_GDS, function(x){x@dataTable@table[, 4] %>% as.numeric} )  %>% as.data.frame
    colnames(tmpdf) <- sapply(all_GDS, function(x){x@header$dataset_id[1]} )
    tmpdf <- tmpdf %>% gather(gds,eval,starts_with('GDS'))
    pdf("/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/EDA/expValues-logOrNonLog1stEd.pdf")
    for( g in unique(tmpdf$gds) ) {
        p <- ggplot(tmpdf %>% filter(gds == g )) + geom_density(aes(x=eval), fill="red", alpha=.5) + labs(title=g)
        print(p)
    }
    dev.off()
}

# Note: some rows like 1320_at or 1405_i_at in GDS4522 are <NA>.
# 
GDS_on_GPL570v <- sapply( all_GDS, function(gds) {gds@header$dataset_id[1]} )

# all GDS in all_GDSv               ~ 2.5 GB
# only GDS with the platform GPL570 ~ 1.2 GB

# Extract all to ESETs
all_ESET <- list()
for (i in seq_along(all_GDS)) {
    all_ESET[[i]] <- GDS2eSet(GDS=all_GDS[[i]], do.log2 = TRUE)
    
    if( any(is.nan(exprs(all_ESET[[i]]))) ){
        # Note: log2(x) is undefined and return NaN if x is negative.
        # What's worse, some datasets looks like already applied log2 transformaiton.
        message("    ", all_GDS[[i]]@header$dataset_id[1], " seems already applied log2(). SKIP")
        all_ESET[[i]] <- GDS2eSet(GDS=all_GDS[[i]], do.log2 = FALSE)
    }
}

function(){
    expl <- sapply(all_ESET, function(x){ exprs(x) } )
    pdf("/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/EDA/EsetVal1stEd.pdf")
    l <- expl[[1]]
    for( l in expl ) {
        print(hist(l))
    }
    dev.off()
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
    all_M[[i]] <- all_M[[i]] %>% na.omit
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
sapply(all_M, function(x) {sum(is.na(x))} )  # MAYBE-LATER What's these nas?
sapply(all_GDS, function(x) {message(paste0(unique(x@dataTable@columns$disease.state),collapse="\t"))} )  # MAYBE-LATER What's these nas?

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
function(){ # Pre-Analysis: What experimental condition each GDS has?
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
}

selectConditionv <- function( m, g = 'GDS5204' ) { # MAYBE-LATER not implemented yet
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

out = list()
i <- k <- 3
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

    dMatrix <- model.matrix( ~ g@dataTable@columns$disease.state )
    
    lmfitted <- limma::lmFit(m, design = dMatrix) # plot(lmfitted$coefficients, pch=18) what's this?
    ebayesed <- eBayes(lmfitted)
    tabled <- topTable(ebayesed, number = 1000)#, p.value = .05)
    tabled
    out[[k]] <- list( lm = lmfitted, table = tabled )
}

tmp_rld <- out[[1]]$table %>% add_rownames('gene') %>% select(gene)
for (k in seq_along(all_M)) {
    options(scipen = 10000)
    m <- all_M[[k]]
    
    this <-
        out[[k]]$table %>% select(P.Value, adj.P.Val) %>%
        add_rownames('gene') %>% mutate_each(funs( r = round(.,5)*100), P.Value, adj.P.Val) %>%
        unite_(col = paste0(attr(m, 'GDS'),'_percent'), c('P.Value', 'adj.P.Val'), sep=" -> ")
    
    tmp_rld <- full_join(tmp_rld, this, by="gene")
}
related_genes_df <- tmp_rld %>% unite(str, starts_with('GDS'), remove = FALSE) %>% filter( nchar(str) > 25 ) %>% select(-str)
#table(sapply(out, function(x){rownames(x$table)}))

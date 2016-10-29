#' Set global variables 'gcl' (global constant list).
#' 
#' @param rootDir Our project (BraindDiseaseCors) directory path. Default is mine.
#' @return Nothing.
setGlobalConstantList <- function(
  rootDir = '/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors'
){
    user <- Sys.info()['user'] 
  if ( user == 'PCUser' ) {
      if ( Sys.info()['sysname'] == 'Darwin' ){
          rootDir <- '/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors'
      } else {
          rootDir <- '/home/PCUser/BrainDiseaseCors'          
      }
  } else {
      rootDir <- rootDir
  }
  message('Root Directory is set to ', rootDir)
  options(max.print=1000) # prevent print flooding
  assign('gcl',
         list(
           affy2uni = paste0(rootDir, '/caprice/MAP/affy2uni.txt'),
           affy2gene = paste0(rootDir, '/caprice/MAP/affy2gene.txt'),
           isPCUser = ifelse(Sys.info()['user'] == 'PCUser', TRUE, FALSE) # to use my cache folder
         ),
         envir = .GlobalEnv
  )
}

#' Set global variables 'gcl' (global constant list).
#' 
#' @param rootDir Our project (BraindDiseaseCors) directory path. Default is mine.
#' @return Nothing.
loadLibraries <- function(
){
  # If you don't have the following libraries,
  # Please download from https://www.bioconductor.org/install/
  library(GEOquery)
  library(limma)   # as.matrix.ExpressionSet is defined in limma
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
}

#' Download all datasets from Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/).
#' 
#' @param GDSnumberv A GDS (Geo Datasets) number vector to be downloaded.
#' @param skipv A skip condition vector.
#' @return A list of downloaded datasets (class 'GDS').
#' @examples
#' download_GDSs(skipv='not-GPL570')
download_GDSs <- function(
  GDSnumberv = c(5204,4879,4859,4838,4758,4532,4522,4523,4477,4414,
                 4358,4231,4218,4154,4136,4135,3834,3502,3459,3345,
                 3129,3128,3113,3110,3069,2978,2941,2821,2795,2613,
                 2191,2190,2154,1962,1917,1912,1835,1816,1815,1813,
                 1726,1253,1096,1085, 910, 909, 833, 810, 707, 596,
                 564, 426, 232, 181),
  skipv = c('not-GPL570', 'no-disease.state', 'no-control',
            'non-binary', 'found-NA')
){
  all_GDS <- list()
  i <- 1
  for ( no in all_GDSv ) {
    if ( gcl$isPCUser ){
      # To avoid downloading the same data again and again.
      # Sometimes downloaded data were corrupted without any warning. Scary
      thisGDS <- getGEO(GEO=paste0("GDS", no), destdir = '/Users/PCUser/Downloads/Rtmp')
    } else {
      thisGDS <- getGEO(GEO=paste0("GDS", no))
    }
    #Sys.sleep(5)
    if( 'not-GPL570' %in% skipv && ! thisGDS@header$platform %in% c('GPL570') ) {
      message('----------Skip GDS ', no, ' because of platform')
      next
    }
    if( ! 'disease.state' %in% colnames(thisGDS@dataTable@columns) ){
        if ('no-disease.state' %in% skipv ){
            message('----------Skip GDS ', no, ' because of no disease.state')
            next
        }
    } else {
        thisGDS@dataTable@columns$disease.state <- gsub(' ','_',thisGDS@dataTable@columns$disease.state)
    }
    if( 'no-control' %in% skipv && (! any(grepl('control', thisGDS@dataTable@columns$disease.state))) ){
      message('----------Skip GDS ', no, ' because of no control')
      next
    }
    if( 'non-binary' %in% skipv && length(unique(thisGDS@dataTable@columns$disease.state)) != 2 ){
      message('----------Skip GDS ', no, ' because of control values')
      next
    }
    if( 'found-NA' %in% skipv && any(is.na(thisGDS@dataTable@table)) ) {
      message('----------Skip GDS ', no, ' because of NA')
      next
    }
    all_GDS[[i]] <- thisGDS
    i <- i + 1
  }
  # all 54 GDS in our 1st proposal    ~ 2.5 GB
  # only GDS with the platform GPL570 ~ 1.2 GB
  return(all_GDS)
}

#' Convert GDS class objects into ExpressionSet class objects.
#' 
#' @param GDSl A list of 'GDS' class objects.
#' @return A list of 'ExpressionSet' class objects.
#' @examples
#' convertGDS2ESET(GDSl)
convertGDS2ESET <- function(
  GDSl
){
  all_ESET <- list()
  for (i in seq_along(GDSl)) {
    
    if ( gcl$isPCUser ){
      # options('download.file.method.GEOquery') == 'auto' by defalt,
      # and libcurls is used for downloading files. Unfortunately, 
      # libcurl does not work on Google Cloud for some reason. So use wget instead
      options('download.file.method.GEOquery'='wget')
      all_ESET[[i]] <- GDS2eSet(GDS=GDSl[[i]], do.log2 = TRUE)
    } else{
      all_ESET[[i]] <- GDS2eSet(GDS=GDSl[[i]], do.log2 = TRUE)
    }
    
    if( any(is.nan(exprs(all_ESET[[i]]))) ){
      # Note: log2(x) is undefined and return NaN if x is negative.
      # What's worse, some datasets looks like already applied log2 transformaiton.
      message("    ", GDSl[[i]]@header$dataset_id[1], " seems already applied log2(). SKIP")
      all_ESET[[i]] <- GDS2eSet(GDS=GDSl[[i]], do.log2 = FALSE)
    }
  }
  # all 54 GDS in our 1st proposal    ~ 1.6 GB
  # only GDS with the platform GPL570 ~ 0.7 GB 
  return(all_ESET)
}

#' Convert GDS class objects into ExpressionSet class objects.
#' 
#' @param ESETl A list of 'ExpressionSet' class objects.
#' @return A list of 'Matrix' class objects.
#' @examples
#' extractMatrixFromEset(ESETl)
extractMatrixFromEset <- function( # FIXME: currently broken.
  ESETl,
  nametype = c('uni', 'gene')[2]
){
  #clip <- pipe("pbcopy", "w"); write.table(unique(affy2uni_df$affy), file=clip, row.names = FALSE); close(clip) # copy to clipboard    
  if ( nametype == 'uni' ){
      renamed_colname <- 'affy'
      affy2uni_df <-
          read.delim(gcl$affy2uni, sep='\t') %>%          # %>% is a pipe. Output of the left is passed to the right
          dplyr::select(-X) %>%                           # avoid confliction between MASS::select() and dplyr::select()
          separate_rows(UniProt.Accession, sep = ';') %>% # separate multiple values like 'Q5TK75; Q6IUU8; Q6V4Z6;'.
          filter( UniProt.Accession != '-' ) %>%          # NB: some Affy ID has no corresponding Uniplot IDs. Strange.
          # This means that our following results are just within some portion of all genes.
          rename( uni = UniProt.Accession ) %>%
          rename_(.dots = setNames('Affy.ID', renamed_colname)) # by using rename_() not rename(), we can use string variables
  } else if ( nametype == 'gene' ) {
      library(hgu133plus2.db) # CAUTION: this library masks dplyr's verbs like select or rename
      library(annotate)
      ls("package:hgu133plus2.db")
      hgu133plus2() # displays when this DB is created
      columns(hgu133plus2.db) # returns all available columns. But too heavy, do not use simultaneously in select()
      all_probes <- keys(hgu133plus2.db, keytype='PROBEID')
      all_entrez <- AnnotationDbi::select(hgu133plus2.db, keys=all_probes, columns="ENTREZID")
      # prbs <- mget(as.character(allGDSl[[1]]@dataTable@table$ID_REF), hgu133plus2ENTREZID)
      # sort(table(unlist(prbs), useNA = 'always'),decreasing = TRUE)[1:10]
      # shows that 12317 probes has no Entrez Gene ID.
      #
      # "In the annotation packages, (by default), we hide probesets that map
      # to more than one gene.  This is because most of the time, you probably
      # don't want anything to do with probes that are not specific."
      # From: https://support.bioconductor.org/p/35647/ 
      # length(hgu133plus2SYMBOL) is the same as nrow(allGDSl[[1]]@dataTable@table)
      mapped_probes <- mappedkeys(hgu133plus2SYMBOL)
      genesym.probeid <- as.data.frame(hgu133aSYMBOL[mapped_probes])
      head(genesym.probeid)
      renamed_colname <- 'affy'
      affy2uni_df <-
          read.delim(gcl$affy2gene, sep='\t') %>%
          dplyr::select(-Species) %>%
          rename( gene = Name ) %>%
          mutate( gene = gsub(' ', '_', gene) )
          rename_(.dots = setNames('AFFYMETRIX_3PRIME_IVT_ID', renamed_colname)) # by using rename_() not rename(), we can use string variables
  } else {
      message('invalid nametype')
      return(1)
  }

  newAttrName <- 'GDS'
  
  all_M <- list()
  i <- 5 # something bad is happening in this dataset, GDS4231
  for (i in seq_along(ESETl)) { 
    tmp <- as.matrix(ESETl[[i]]) # NB: call library(limma); and library(affy) before
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
    attr(all_M[[i]], newAttrName) <- ESETl[[i]]@experimentData@other$dataset_id[1]
    
    #all_DF[[i]] <- as.data.frame(ESETl[[i]]) 
    #
    # Applying as.data.frame is not a good idea.
    # It's because R doesn't allow to use rownames starting with numbers.
    # '1007_s_at' or '121_at' is automatically converted to 'X1007_s_at' or 'X121_at'.
    #
    # It induces some errors
    # '[E]specially when you've been looking at your screen for 20 straight hours' :)
    #
    # You can check by
    # m <- as.matrix(    ESETl[[1]])  # NB: call library(limma); and library(affy) before
    # d <- as.data.frame(ESETl[[1]])
    # data.frame( matrix_name = rownames(m)[1:5], df_name = colnames(d)[1:5] )
    # 
    # I think this is the reason TA uses as.matrix in place of as.data.frame
  }
  # all 54 GDS in our 1st proposal  ~ 600 MB
  # all_M with the platform GPL570  ~ 300 MB
  return(all_M)
}


#' Not implemented yet.
#' 
#' @param m 
#' @param g
#' @return 
#' @examples
#' 
selectConditionv <- function(
  m, g = 'GDS5204'
){ # MAYBE-LATER not implemented yet
  condition <- list()
  if ( g %in% paste0('GDS', c(1917, 1962, 2154, 2795, 2821,
                              3502, 4135, 4136, 4218, 4231,
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

#' Apply t-test
#' 
#' @param Ml A list of 'Matrix'. [ Genes x Samples ]  
#' @param GDSl A list of 'GDS' class objects. length(GDSl) must be the same as length(Ml).
#' @return A list containing t-test results for each dataset.
#' @examples applyTtestToGeneExpressionMatrices(Ml, GDSl)
applyTtestToGeneExpressionMatrices <- function(
  Ml,
  GDSl
){
  out = list()
  i <- k <- 3
  for (k in seq_along(Ml)) {
    
    m <- Ml[[k]]
    g <- GDSl[[k]]
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
  return(out)
}

# TODO convert to a function

# tmp_rld <- out[[1]]$table %>% add_rownames('gene') %>% select(gene)
# for (k in seq_along(all_M)) {
#     options(scipen = 10000)
#     m <- all_M[[k]]
#     
#     this <-
#         out[[k]]$table %>% select(P.Value, adj.P.Val) %>%
#         add_rownames('gene') %>% mutate_each(funs( r = round(.,5)*100), P.Value, adj.P.Val) %>%
#         unite_(col = paste0(attr(m, 'GDS'),'_percent'), c('P.Value', 'adj.P.Val'), sep=" -> ")
#     
#     tmp_rld <- full_join(tmp_rld, this, by="gene")
# }
# related_genes_df <- tmp_rld %>% unite(str, starts_with('GDS'), remove = FALSE) %>% filter( nchar(str) > 25 ) %>% select(-str)
# #table(sapply(out, function(x){rownames(x$table)}))

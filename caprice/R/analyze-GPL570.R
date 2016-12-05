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
      if ( user == 'admin' ){
          rootDir <- '/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree'
      } else if (user == 'ianjohnson' ){
          rootDir <- '/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/' # IS THIS RIGHT?
      } else {
          message('unknown user')
          rootDir <- rootDir
      }
  }
  message('Root Directory is set to ', rootDir)
  options(max.print=1000) # prevent print flooding
  assign('gcl',
         list(
           affy2uni = paste0(rootDir, '/caprice/MAP/affy2uni.txt'),
           isPCUser = ifelse(Sys.info()['user'] == 'PCUser', TRUE, FALSE), # to use my cache folder
           rootDir = rootDir,
           user = user
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
  if( gcl$user == 'ianjohnson' ){
      message('early return')
      return()
  }
    
  # If you don't have the following libraries,
  # Please download from https://www.bioconductor.org/install/
  library(GEOquery)
  library(limma)   # as.matrix.ExpressionSet is defined in limma
  library(annotate) # BiocInstaller::biocLite("annotate")
  library(hgu133plus2.db)

  library(tibble)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
}

#' Download all datasets from Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/).
#' Use local caches if possible.
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
            'non-binary', 'found-NA', 'blacklist'),
  suppress=TRUE
){
  all_GDS <- list()
  i <- 1
  for ( no in GDSnumberv ) {
    if ( gcl$isPCUser ){
      # To avoid downloading the same data again and again.
      # Sometimes downloaded data were corrupted without any warning. Scary
      options('download.file.method.GEOquery'='wget')
      if (suppress) {
        thisGDS <- suppressMessages(
          getGEO(GEO=paste0("GDS", no), destdir = '/Users/PCUser/Downloads/Rtmp')
        )
      } else{
        thisGDS <- getGEO(GEO=paste0("GDS", no), destdir = '/Users/PCUser/Downloads/Rtmp')
      }
    } else {
      if (Sys.info()["user"] == "admin" ) {
        thisGDS <- getGEO(GEO=paste0("GDS", no), destdir = '/Users/admin/Downloads/Rtmp') 
      } else {
        thisGDS <- getGEO(GEO=paste0("GDS", no), destdir = '/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/data')
      }
    }
    #Sys.sleep(5)
    if( 'not-GPL570' %in% skipv && ! thisGDS@header$platform %in% c('GPL570') ) {
      message('----------Skip GDS ', no, ' because of platform')
      next
    }
    gds_num <- thisGDS@header$dataset_id[1]
    
    removeBL <- 'blacklist' %in% skipv
    
    if ( removeBL && gds_num == 'GDS4838') { message("GDS4838 is only various tumors");    next;  }
    if ( removeBL && gds_num == 'GDS4532') { message("GDS4532 is not disease study");      next;  }
    if ( removeBL && gds_num == 'GDS4477') { message("GDS4477 is only diseased tumor");    next;  }
    if ( removeBL && gds_num == 'GDS2154') { message("GDS2154 is not brain disease");    next;  }
    if ( removeBL && gds_num %in% c('GDS4231', 'GDS3502')) { message("GDS4231/3502: suspicious"); next;  }
    
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
  GDSl, suppress=TRUE
){
  all_ESET <- list()
  for (i in seq_along(GDSl)) {
    
    if ( gcl$isPCUser ){
      # options('download.file.method.GEOquery') == 'auto' by defalt,
      # and libcurls is used for downloading files. Unfortunately, 
      # libcurl does not work on Google Cloud for some reason. So use wget instead
      options('download.file.method.GEOquery'='wget')
      if (suppress) {
        all_ESET[[i]] <- suppressMessages(GDS2eSet(GDS=GDSl[[i]], do.log2 = TRUE))
      } else {
        all_ESET[[i]] <- GDS2eSet(GDS=GDSl[[i]], do.log2 = TRUE)
      }
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
  nametype = c('uni', 'gene', 'none')[3]
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
          dplyr::rename( uni = UniProt.Accession ) %>%
          dplyr::rename_(.dots = setNames('Affy.ID', renamed_colname)) # by using rename_() not rename(), we can use string variables
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
          dplyr::rename( gene = Name ) %>%
          mutate( gene = gsub(' ', '_', gene) )
          dplyr::rename_(.dots = setNames('AFFYMETRIX_3PRIME_IVT_ID', renamed_colname)) # by using rename_() not rename(), we can use string variables
  } else if ( nametype == 'none' ) {
    # do nothing
  } else {
      message('invalid nametype')
      return(1)
  }

  newAttrName <- 'GDS'
  
  all_M <- list()
  i <- 5 # something bad is happening in this dataset, GDS4231
  for (i in seq_along(ESETl)) { 
    if ( nametype != 'none' ){
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
    } else {
      all_M[[i]] <- as.matrix( ESETl[[i]] )  # e.g. [ 65107 genes x 41 SaMples ]
    }
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

addAgeColumn <- function(
  GDSl
){

}

checkExperimentalCondition <- function(
  GDSl
){
  gdsv <- sapply(GDSl, function(x) {x@header$dataset_id[1]})
  condv <- sapply(GDSl, function(x){ colnames(x@dataTable@columns) } ) %>% unlist %>% unique
  condv <- condv[! condv %in% c('sample', 'description')]
  out <- matrix(NA, nrow=length(GDSl), ncol=length(condv), dimnames = list(gdsv, condv))
  i <- k <- 1
  for( i in seq_along(GDSl) ){
    for( k in seq_along(condv) ){
      if ( condv[k] %in% colnames(GDSl[[i]]@dataTable@columns) ){
        out[i, k] <- length(unique(GDSl[[i]]@dataTable@columns[, condv[k]]))
      }
    }
  }
  return( out )
}

convertDiseaseStateIntoBinary <- function(n = 'GDS5204', gds){
    if ( n == 'GDS5204' ) { # aging study
        # age[4 category], gender is the experimental condition
        # gender_ <- gds@dataTable@columns$gender
        dz_ <- ifelse(gds@dataTable@columns$age%in%c('young (<40yr)', 'middle aged (40-70yr)'), FALSE, TRUE)
        #} else if (n == 'GDS4838') { dz_ <- gds@dataTable@columns$disease.state != 'CNS_primitive_neuroectodermal_tumors' # FIXME: read paper
    } else if (n == 'GDS4523') { dz_ <- gds@dataTable@columns$disease.state                       # ignore age and gender
    } else if (n == 'GDS4522') { dz_ <- gds@dataTable@columns$disease.state                       # ignore age and gender
    } else if (n == 'GDS4358') { dz_ <- gds@dataTable@columns$disease.state != 'control'
    #} else if (n == 'GDS4231') { dz_ <- gds@dataTable@columns$disease.state                       # FIXME: consider therapy # suspicious
    } else if (n == 'GDS4218') { dz_ <- gds@dataTable@columns$disease.state != 'healthy_control'  # MAYBE-LATER only 6 samples
    } else if (n == 'GDS4136') { dz_ <- gds@dataTable@columns$disease.state != 'control'          # Alzheimer
    } else if (n == 'GDS4135') { dz_ <- gds@dataTable@columns$disease.state != 'Braak_stage_I-II' # Alzheimer (stage)
    #} else if (n == 'GDS3502') { dz_ <- gds@dataTable@columns$disease.state != 'control'          # FIXME: remove bipolar_disorder # suspicious 
    } else if (n == 'GDS2821') { dz_ <- gds@dataTable@columns$disease.state != 'control'          # Parkinson
    } else if (n == 'GDS2795') { dz_ <- gds@dataTable@columns$disease.state != 'normal'           
    #} else if (n == 'GDS2154') { dz_ <- gds@dataTable@columns$disease.state != 'healthy'
    } else if (n == 'GDS1962') { dz_ <- gds@dataTable@columns$disease.state != 'non-tumor'
    } else if (n == 'GDS1917') { dz_ <- gds@dataTable@columns$disease.state != 'control'
    } else {
        message("GDS: ",n )
        stop("not implemented")
    }
    return( dz_ )
}

getGenderMetaInfo <- function(targetGDS='GDS4358'){
    # see get_gender.R
    library(GEOmetadb) # not in loadLibraries()
    
    # getSQLiteFile()
    # If the above function getSQLiteFile() fails, just download from
    # https://dl.dropboxusercontent.com/u/51653511/GEOmetadb.sqlite.gz
    
    if ( Sys.info()['user'] == 'PCUser' ) {
        metadb_path <- '/Users/PCUser/Downloads/Rtmp/GEOmetadb.sqlite'
    } else {
        metadb_path <- '/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree/zangcc/GEOmetadb.sqlite'
    }
    con <- dbConnect(SQLite(), metadb_path)
    geo_tables <- dbListTables(con)
    # geo_tables
    # dbListFields(conn = con, name = 'gse_gsm')
    # sapply(geo_tables, function(x) dbListFields(con,x))
    
    # dbGetQuery(con, "select * from gds limit 3")
    # You can use dbGetQuery here, but you can also access to GEOmetadb.sqlite file by dplyr
    
    db <- src_sqlite(metadb_path)
    GSEv <- tbl(db,'gds') %>% filter( gds == targetGDS ) %>%
        dplyr::select(gse) %>% collect() %>% c %>% unlist
    
    #Do some formating like adding quote for characters to create querries
    partialQuery <- paste0("'", paste(GSEv, collapse="', '"), "'" )
    
    gsmdf <- dbGetQuery(con, paste("select * from gse_gsm where gse IN(", partialQuery, ')'))
    tmp <- tbl(db, 'gsm') %>% filter( gsm %in% gsmdf$gsm ) %>%
        dplyr::select(gsm, characteristics_ch1) %>% collect
    
    genderdf <-
        tmp %>% filter( 
            (!grepl('[sS]tage|[dD]osage|[lL]ineage|[pP]assage',characteristics_ch1) &
                 grepl('([sS]ex|[gG]ender)', characteristics_ch1))) %>%
        dplyr::rename( gender = characteristics_ch1) %>%
        mutate( gender = str_extract(string=gender,
                                     pattern = '([sS]ex|[gG]ender): ([mM]|[mM]ale|[fF]|[fF]emale)') ) %>%
        mutate( gender = gsub('.* ','', gender) %>% toupper )

    genderdf %>% mutate( gsm = as.factor(gsm) ) %>% dplyr::rename( sample = gsm )
}

#' Not implemented yet.
#' 
#' @param m 
#' @param g
#' @return 
#' @examples
#' 
buildDesignMatrix <- function(
  gds, type = c('binary', 'gender')[1]
){ # MAYBE-LATER not implemented yet
  
  n <- gds@header$dataset_id[1]
  if ( type == 'binary' ){
    dz_ <- convertDiseaseStateIntoBinary(n, gds)
    dMatrix <- model.matrix( ~ dz_ )
  } else if (type == 'gender') {
      # GDS5204, GDS4522, GDS4523, GDS2821 has gender columns in dataTable@columns
      # GDS4358, GDS4136 has gender information in meta data
      dz_ <- convertDiseaseStateIntoBinary(n, gds)
      if ( n %in% c('GDS5204', 'GDS4522', 'GDS4523', 'GDS2821')) {
          # these four datasets denotes females as "female"
          isFemale_ <- grepl('female', gds@dataTable@columns$gender)
          dMatrix <- model.matrix( ~ dz_ + isFemale_ )
      } else if (n %in% c('GDS4358', 'GDS4136')) {
          genderdf <-
              left_join(gds@dataTable@columns,
                        getGenderMetaInfo(n),
                        by = 'sample')
          isFemale_ <- grepl('F', genderdf$gender) # both store as 'F'
          dMatrix <- model.matrix( ~ dz_ + isFemale_ )
      } else {
          dMatrix <- model.matrix( ~ dz_ )
      }
      
  }
  
  return(dMatrix)
}

#' Apply t-test
#' 
#' @param Ml A list of 'Matrix'. [ Genes x Samples ]  
#' @param GDSl A list of 'GDS' class objects. length(GDSl) must be the same as length(Ml).
#' @return A list containing t-test results for each dataset.
#' @examples applyTtestToGeneExpressionMatrices(Ml, GDSl)
applyTtestToGeneExpressionMatrices <- function(
  Ml,
  GDSl,
  nTopGene = 1000,
  p_threshold = NULL,
  type = c('fixed_n', 'p', 'otherwise')[1],
  method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[7],
  design = c('binary', 'gender')[1],
  removeGenderEffect = FALSE
){
  out = list()
  i <- k <- 3
  message("build ", design, " design matrices.")
  for (k in seq_along(Ml)) {
    
    m <- Ml[[k]]
    g <- GDSl[[k]]
    gdsno <- g@header$dataset_id[1]
    
    if ( gdsno == 'GDS4358' ){
        if ( removeGenderEffect ){
            message("removeGenderEffect from GDS4358")
            genderdf <-
                left_join(g@dataTable@columns,
                          getGenderMetaInfo(gdsno),
                          by = 'sample')
            g@dataTable@columns <- genderdf %>% filter( gender != 'F' )
            m <- m[, g@dataTable@columns$sample]
        }
        # It seems that three different regions are sampled from each patient.
        # Here we chose frontal cortex due to compatibility with other datasets
        g@dataTable@columns <- g@dataTable@columns %>% filter( tissue == 'Frontal cortex' )
        m <- m[, g@dataTable@columns$sample]
        message("frontal cortex is selected in GDS4358")
    }
    
    dMatrix <- buildDesignMatrix(g, type = design) # FIXME: construct design matrices for each dataset
    
    message(attr(m, 'GDS'), appendLF = FALSE)
    if ( is.null(dMatrix) ){
      message('  SKIP unimplemented design matrix')
      next  
    }else {
      message()
    }
    
    
    # as.numeric(gsub('[^0-9]','',GDSl[[1]]@dataTable@columns$age))
    # model.matrix( ~ GDSl[[1]]@dataTable@columns$disease.state + as.numeric(gsub('[^0-9]','',GDSl[[1]]@dataTable@columns$age)) )
    
    lmfitted <- suppressMessages(limma::lmFit(m, design = dMatrix)) # plot(lmfitted$coefficients, pch=18) what's this?
    ebayesed <- eBayes(lmfitted)
    if ( type == "fixed_n" ) { # FIXME 400 and more
        #tabled <- topTable(ebayesed, p.value = p_threshold, adjust.method = method)
        #if ( nrow(tabled) < 1 ) {
        #    message(gdsno, " too few degs ... return ",nTopGene," tops instead")
            tabled <- topTable(ebayesed, number = nTopGene, adjust.method = method)
        #}
    } else if (type == 'p') {
        tabled <- topTable(ebayesed, p.value = p_threshold, adjust.method = method)
    } else {
        stop("not implemented")
    }
    if( nrow(tabled) == 0 ){
        tabled <- data.frame(matrix(0, nrow=1, ncol=6))
        colnames(tabled) <- c("logFC", "AveExpr", "t",
                              "P.Value", "adj.P.Val", "B")
    }
    out[[k]] <- list( lm = lmfitted, table = tabled, dm = dMatrix, gds = gdsno )
  }
  return(out)
}

#' Select out significant genes
#' 
#' @param Ml A list of 'Matrix'. [ Genes x Samples ]  
#' @param GDSl A list of 'GDS' class objects. length(GDSl) must be the same as length(Ml).
#' @return A list containing t-test results for each dataset.
#' @examples applyTtestToGeneExpressionMatrices(Ml, GDSl)
pickSignificantGenes <- function(
  out,
  Ml
){
  tmp_rld <- out[[1]]$table %>% rownames_to_column('probe') %>% dplyr::select(probe)
  for (k in seq_along(Ml)) {
    options(scipen = 10000) # for avoiding scientific digit notation
    m <- Ml[[k]]
    no <- gsub('GDS', '', attr(m, 'GDS'))
    this <-
      out[[k]]$table %>% dplyr::select(P.Value, adj.P.Val) %>%
      rownames_to_column('probe') %>% mutate( P.Value = round(P.Value,5)) %>%
      mutate( adj.P.Val = round(adj.P.Val,5)) %>%
      dplyr::rename( pval = P.Value, adjpval = adj.P.Val) %>%
      unite_(col = paste0('adjust'), c('pval', 'adjpval'), sep=" -> ", remove = FALSE) %>%
      set_names( c('probe', paste0(names(.)[names(.)!='probe'], '_', no)))
    
    tmp_rld <- full_join(tmp_rld, this, by="probe")
  }
  related_genes_df <- tmp_rld %>% unite(str, starts_with('adjust'), remove = FALSE) %>%
    filter( nchar(str) > 25 ) %>% dplyr::select(-str, -starts_with('adjust'))
  return(related_genes_df)
}


build_deg_matrix <- function(
    topped, type = c('t', 'F')[1]
){
    
    type_up <- 'up'
    type_lo <- 'lo'
    
    degdf <- 
        lapply(topped, function(x) {
        tmpgds <- x$gds;
            if( type == 't' ) {
                x$table %>% rownames_to_column("probe")  %>%
                    mutate( gds = tmpgds, type = ifelse(t > 0, type_up, type_lo) )
            } else {
                x$table %>% rownames_to_column("probe")  %>%
                    mutate( gds = tmpgds, type = type_lo ) # type doesn't matter
            }
        } ) %>% rbindlist
    
    gdsnov <- sapply(topped, function(x) { x$gds })
    
    k <- 1
    i <- 2
    
    deg_matrix <- matrix( NA, ncol = length(topped), nrow = length(topped) )
    
    rownames(deg_matrix) <- gdsnov
    colnames(deg_matrix) <- gdsnov
    
    for( i in seq_along(topped) ){
        for( k in i:length(topped) ){
            # AveExpr == 0 is a placeholder used in applyTtestToGeneExpressionMatrices
            if( topped[[k]]$table$AveExpr[1] == 0 ||
                topped[[i]]$table$AveExpr[1] == 0 ) {
                    deg_matrix[i,k] <- 0
                    deg_matrix[k,i] <- 0
                next
            }
            
            gds_i <- gdsnov[[i]]
            gds_k <- gdsnov[[k]]
            
            deg_i_up <- (degdf %>% filter( gds == gds_i & type == type_up ))$probe
            deg_i_lo <- (degdf %>% filter( gds == gds_i & type == type_lo ))$probe
            deg_k_up <- (degdf %>% filter( gds == gds_k & type == type_up ))$probe
            deg_k_lo <- (degdf %>% filter( gds == gds_k & type == type_lo ))$probe
            
            if( i == k )
                next
            deg_matrix[i,k] <- length(intersect(deg_i_up, deg_k_up))
            deg_matrix[k,i] <- length(intersect(deg_i_lo, deg_k_lo))
        }
    }
    attr(deg_matrix, 'df') <- degdf
    deg_matrix
}

#' perform Over Representation Analysis (one of Pathway Enrichment Analysis)
#' 
#' @param relatedGenev Usually a vector of rownames(topTable())
#' @param gene2pathwaydf parsed Uniprot2Reactome.txt data frame
#' @return a data frame with 2 columns (pathwayName, p-value)
perform_OverRepresentationAnalysis <- function( relatedGenev, gene2pathwaydf, n_min = 2, n_max = 500) {
    
    relatedGenev <- unique(relatedGenev)
    
    # I did ORA in a different way from TA's script.
    #
    # First, I defined two matrices as follows:
    #   M_PATH_x_ALLGENES  : the [i, j] entry is 1, if the j-th gene are in that i-th pathway. Otherwise 0.
    #   M_ALLGENES_RELATED : the [i, j] entry is 1, if the j-th related gene is the same as the i-th gene.
    #
    # Then, the dot product of M_PATH_x_ALLGENES and M_ALLGENES_RELATED are calculated into M_PATH_x_RELATED.
    # M_PATH_x_RELATED[i, j] is 1 if the j-th related genes are involved in that i-th pathway.
    #
    # All the varialbes defined below are for that matrix computations. So feel free to skip them.
    #
    # This approach is much faster than using apply() like TA described.
    # And, I use Matrix::sparseMatrix data structure for the sake of memory efficiency.
    
    # Getting background genes (all the unique uniprot genes in this file) increases
    # Sensitivity (True Positive Rate) and Specificity (True Negative Rate) for Pathway Enrichment Analysis.
    # It's possible to perform PEA only with interested genes, but it decreases TPR and TNR.

    allGenev <- sort(unique(gene2pathwaydf$uni))
    # 10467 for both UniProt2Reactome.txt and UniProt2Reactome_All_Levels.txt
    
    # The Uniprot2Reactome.txt TA gave to us contains non-homo-sapies uniprots as well. Should we?
    
    candidatePathwayv <-
        gene2pathwaydf %>%
        group_by(pathway) %>%
        summarize( n = n() ) %>%            # calculate number of genes (rows) for each reactome pathway
        filter( n_min < n & n < n_max ) %>% # only keep the pathways with moderate number of genes
        .$pathway %>% sort                  # '.' (dot) stores the previous output (you can regard it as stdin).
    
    head(candidatePathwayv) 
    # => If success, we get 1136 pathways for UniProt2Reactome.txt,
    #                       1516 pathways for UniProt2Reactome_All_Levels.txt
    
    pathwaydf <- gene2pathwaydf %>% filter( pathway %in% candidatePathwayv ) 
    
    path2r        <- c( 1:length(pathwaydf$pathway %>% unique %>% sort) ) # r is row index
    names(path2r) <-             pathwaydf$pathway %>% unique %>% sort 
    
    gene2c        <- c( 1:length(pathwaydf$uni %>% unique %>% sort) )     # c is column index
    names(gene2c) <-             pathwaydf$uni %>% unique %>% sort
    
#   path2r[pathwaydf$pathway]
#    gene2c[pathwaydf$uni]
    
    # [ numberOfPathways x numberOfGenes ] = [ 2223 x 85875 ],
    # But it's really sparse since we filtered pathways.
    # (In each row, there are at most 500 entries)
    M_PATH_x_ALLGENES <-
        Matrix::sparseMatrix( dims = c(length(candidatePathwayv), length(allGenev)),
                              i = path2r[pathwaydf$pathway],
                              j = gene2c[pathwaydf$uni], 
                              x = 1, # values for non-zero entries. Use 1 (multiplicative identity)
                              symmetric = FALSE, triangular = FALSE,
                              index1 = TRUE # row and column indexes start from 1, not 0
        )
    
    # some of the related genes are not in gene2c for two reasons:
    # 1. because we filtered n_min < n < n_max
    # 2. because UniProt2Reactome.txt does not contain entries for all UniProt Accessions
    enrichedGenev <- gene2c[ relatedGenev[ ! is.na(gene2c[relatedGenev]) ] ]
    message("Among ", length(relatedGenev), " Uniprots, ",
            length(enrichedGenev), " has entries in Reactome DB.")
    message("Over-representation Analysis -- Among ", length(unique(gene2pathwaydf$pathway)), " Pathways, ",
            length(unique(pathwaydf$pathway)), " pathways will be searched.")
    
    M_ALLGENES_x_RELATED <-
        Matrix::sparseMatrix( dims = c(length(allGenev), length(enrichedGenev)),
                              i = enrichedGenev,
                              j = 1:length(enrichedGenev), 
                              x = 1, # values for non-zero entries. Use 1 (multiplicative identity)
                              symmetric = FALSE, triangular = FALSE,
                              index1 = TRUE # row and column indexes start from 1, not 0
        )
    
    #  [ 1516 x 166 ]     [ 1516 x 10467 ]         [ 10467 x 166  ]
    M_PATH_x_RELATED <- M_PATH_x_ALLGENES %*% M_ALLGENES_x_RELATED
    colnames(M_PATH_x_RELATED) <- names(enrichedGenev)
    
    message("M_PATH_x_RELATED (", nrow(M_PATH_x_RELATED), " x ", ncol(M_PATH_x_RELATED) ,")",
            " have been built successfully.")
    
    # ng_ : Number of Genes
    ng_all       <- length(allGenev)                                       # 10467
    # ng_related   <- length(relatedGenev)                                 #  454
    ng_related   <- length(enrichedGenev)                                  #  166
    ng_enrichedv <- Matrix::rowSums(M_PATH_x_RELATED,  sparseResult=FALSE) #  1516 elements
    ng_pathv     <- Matrix::rowSums(M_PATH_x_ALLGENES, sparseResult=FALSE) #  1516 elements
    
    # I decided to use length(enrichedGenev) as ng_related
    # since it's the number of columns of M_PATH_x_RELATED and more convincing

    # Let 
    #          ng_all as the number of     all genes (allGenev)
    #      ng_related as the number of related genes (relatedGenev)
    #     ng_pathv[i] as the number of genes in the i-th pathway (candidatePathwayv[i])
    # ng_enrichedv[i] as the number of genes related to the disease and in the i-th pathway (ng_enrichedv[i])
    #
    # Then construct the following 2 x 2 table for pathway[i]:
    #
    #                                  related                                       Non-related
    #     in-Pathway[i]                   ng_enrichedv[i]                             ( ng_pathv[i] - ng_enrichedv[i]  )
    # not-in-pathway[i]    ( ng_related - ng_enrichedv[i] )    ( ng_all - (ng_related + ng_pathv[i] - ng_enrichedv[i]) )
    #
    # For details, see: p.207 of https://www.ncbi.nlm.nih.gov/pubmed/23192548
    # In the above book,
    #                     related       Non-related            total
    #     in-Pathway[i]        k              m - k                m
    # not-in-pathway[i]    n - k     N - (n + m - k)           N - m
    #             total    n         N -  n                    N
    #

    out <-
        data.frame( pathway = candidatePathwayv, n_path = ng_pathv, n_enriched = ng_enrichedv ) %>%
        mutate( pval = NA, fdr = NA, gene = NA)
    
    i <- 1  # This i is for test purpose. Since it is outside of for loop, no effect in production 
    for( i in 1:nrow(out) ){
        cont_table <- # contingency table for each pathway
            matrix(
                data =
                    c(              ng_enrichedv[i],                           ng_pathv[i] - ng_enrichedv[i],
                       ng_related - ng_enrichedv[i],    ng_all - (ng_related + ng_pathv[i] - ng_enrichedv[i]) 
                    ),
                nrow = 2, ncol = 2
            )
        # 3. For each pathway, use the fisher.test() to assess the enrichment of DiffExpGenes in the pathway.
        # We need is one-side test and set the altertive hypothesis as "greater". 
        out[i, 'pval'] <- fisher.test(cont_table, alternative = c('two.sided', 'greater', 'less')[2])$p.value
        out[i, 'gene'] <- paste0(colnames(M_PATH_x_RELATED)[M_PATH_x_RELATED[i, ] == 1], collapse='|')
        # validated by
        # out %>% mutate( nn = str_count(gene, '\\|') + 1) %>% filter( gene != '') %>% filter( nn != n_enriched)
    }
    out[, 'fdr']  <- p.adjust(out[, 'pval'], method="fdr");
    out
}

do_pea <- function( gene2pathwaydf, probev, p_threshold = .1, nopath=FALSE, n_min=2, n_max=1000 ){
    probe2uni <- suppressMessages(
                    AnnotationDbi::select(hgu133plus2.db, keys=probev,
                                           columns="UNIPROT") %>% filter( ! is.na(UNIPROT) )
                                )
    relatedGenev <- unique(probe2uni$UNIPROT)
    message("Among ",length(probev)," given probes, ",
            length(unique(probe2uni$PROBEID)), " probes have ",
            length(relatedGenev), " Uniprot Accessions in hgu133plus2.db database.")
    
    oraed <- perform_OverRepresentationAnalysis( relatedGenev, gene2pathwaydf,
                                                 n_min = n_min, n_max = n_max)
    if (nopath)
        oraed %>% filter( pval <= p_threshold ) %>% dplyr::select(-n_path)
    else
        oraed %>% filter( pval <= p_threshold )
}


# Note: this is for Limma. SAM should be handled separately
do_PathwayEnrichmentAnalysis <- function(gene2pathwaydf, topped, gdsv, gds = 'GDS5204', p_threshold = 1, nopath=FALSE){
    # p_threshold == 1 returns all pathways.
    # gdsv <- sapply(topped, function(x) x$gds )
    probev <- rownames( topped[[ which(gdsv == gds) ]]$table )  
    peaed <- do_pea(gene2pathwaydf, probev, p_threshold = 1, nopath=nopath)
    peaed
}

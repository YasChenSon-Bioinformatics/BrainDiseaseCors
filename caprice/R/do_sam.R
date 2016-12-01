source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')
library('limma')


source('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis/functions.R')
setwd('/Users/ianjohnson/Desktop/Columbia/bioinformatics/project/sam_analysis')

datasets_num = list(4523, 4522, 4358, 4218, 4136,
                    2821, 1917, 5204, 4135, 2795,
                    1962)
datasets = download_GDS(datasets_num)



# COMPUTE SAM
#
# SUPER HEAVY ( ~ 17 GB). Use a generated csv file below instead
#
bool_to_num = function(x){ if (x==TRUE) 1 else 2}
all_sam = list()
i = 1
for (df in datasets) {
    message(datasets_num[[i]], ' -----------------------------------------------')
    
    # if (i < 8) { i = i + 1; next; } # uncomment to fast forward
    #if (i == 4) {
    if (datasets_num[[i]] %in% c(4522, 4218, 2821, 1962) ) {
        # more robust for future change
        
        #i = i + 1; next
        # table(colSums(is.na( df )) == nrow(df) )
        # FALSE  TRUE 
        # 33398 21277 
        # Almost 2/5 of probes are all NA.
        # Remove the probe columns in order to avoid the following SAM() errors:
        #  'UWXYZ rows with more than 50 % entries missing'
        #
        # sapply(datasets,
        #        function(df){ nrow_allNA = sum(colSums(is.na( df )) == nrow(df)); nrow_allNA })
        # returns 0  2728     0 21277     0   357     0     0     0     0    62
        #
        isAllNAcolv <- colSums(is.na( df )) == nrow(df)
        df <- df[ , ! isAllNAcolv ]
        message("rows with all NA values are removed.")
        # Still, 4522 causes errors because of almost-all NA rows. 
    } # error in 4th dataset...
    
    # Possible levels used for control are:
    # 
    # "control"
    # "normal" (2795)
    # "young" (aging study 5204)
    # "Braak stage I-II" (alz 4135)
    # "non-tumor" (scf/angiogenesis 1962)
    
    control_strings = c("control", "normal", "Braak stage I-II", "non-tumor")
    dz_ctrl_boolv = grepl(paste(control_strings, collapse = "|"), df$disease.state)
    
    # Aging study has no disease.state column - use "age"
    # Note - grepl doesn't seem to like special characters
    if (datasets_num[[i]] == 5204) {
        age_control_strings = c("young", "middle")
        dz_ctrl_boolv = grepl(paste(age_control_strings, collapse = "|"), df$age)
    }
    
    y = apply(as.data.frame(dz_ctrl_boolv), 1, bool_to_num)
    
    # i == 4 (GDS4218) has remained a "genotype/variation" column and caused clash
    # checked by rownames(df_t)[!grepl("[0-9]",rownames(df_t))]
    drops <- c("sample", "age", "gender", "tissue", "genotype.variation",
               "genotype/variation",
               "development.stage", "agent", "other", "cell.type", "disease.state",
               "description", "individual")
    df = df[, !(names(df) %in% drops)]
    df_t = t(df)
    
    message(  paste0(rownames(df_t)[!grepl("[0-9]",rownames(df_t))], collapse=" ") )
    
    samfit <- SAM(df_t, y, resp.type="Two class unpaired", nperms=1000, fdr.output = 0.01,
                  geneid = rownames(df_t) )
    all_sam[[i]] = samfit
    i = i + 1
}
#########################################################################################




# PLOT SAM ##############################################################################
par(mfrow=c(3,3))
par(cex.axis=1.5, cex.lab=1.75, cex.main=1.5, cex.sub=1.5)
for (sam in all_sam) {
    plot(sam)
}


# Collect DEGs ##############################################################################
i = 1
deg_up = list()
deg_lo = list()
for (sam in all_sam) {
    message(i, typeof(all_sam[[i]]$siggenes.table$genes.up), " ", typeof(all_sam[[i]]$siggenes.table$genes.lo))
    #if (i == 4) { i = i + 1; next } # error in 4th dataset...
    
    # Take the Gene Names (e.g. 6473)
    deg_up[[i]] = all_sam[[i]]$siggenes.table$genes.up[, 'Gene Name']
    deg_lo[[i]] = all_sam[[i]]$siggenes.table$genes.up[, 'Gene Name']
    
    i = i + 1
}

# Save CSVs
for (i in 1:length(deg_up)) {
    if (Sys.info()['user'] == 'PCUser'){
        library(dplyr)
        deg <- list()
        
        type_up <- 'up'
        type_lo <- 'lo'
        
        uped <- rbindlist(
            sapply(all_sam, function(sam) sam$siggenes.table$genes.up  %>% as.data.frame )[
                which(sapply(all_sam, function(sam) ! is.null(sam$siggenes.table$genes.up)))
            ])
        uped$gds = rep(unlist(datasets_num), sapply(all_sam, function(x) x$siggenes.table$ngenes.up))
        uped$type = type_up
        
        loed <- rbindlist(
            sapply(all_sam, function(sam) sam$siggenes.table$genes.lo  %>% as.data.frame )[
                which(sapply(all_sam, function(sam) ! is.null(sam$siggenes.table$genes.lo)))
                ])
        loed$gds = rep(unlist(datasets_num), sapply(all_sam, function(x) x$siggenes.table$ngenes.lo))
        loed$type = type_lo
        
        degdf <- bind_rows(uped, loed)
        colnames(degdf) <- c("geneid", 'probe', 'score', 'numerator',
                             'denominator', 'fold-change', 'qval', 'gds', 'type')
        write.csv(degdf, file = "caprice/RESULT/11GDS_sam_results.csv", quote=FALSE, row.names=FALSE)
        
        degdf %>% group_by(probe, type) %>% mutate( n_dataset = n() ) %>% filter( n_dataset > 1 ) %>% summarize( gds = collapse() )
        
        degdf %>% group_by(gds, type) %>% summarize(n = n() )
        
        k <- 1
        i <- 2
        
        deg_matrix <- list()
        
        deg_matrix <- matrix( NA, ncol = length(datasets), nrow = length(datasets) )

        rownames(deg_matrix) <- paste0('GDS', datasets_num) 
        colnames(deg_matrix) <- paste0('GDS', datasets_num)
    
        for( i in seq_along(datasets) ){
            for( k in i:length(datasets) ){
                
                gds_i <- datasets_num[[i]]
                gds_k <- datasets_num[[k]]
                
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
    
        #         GDS4523 GDS4522 GDS4358 GDS4218 GDS4136 GDS2821 GDS1917 GDS5204 GDS4135 GDS2795 GDS1962
        # GDS4523      NA       0       0       0       0       0       0       0       0       0       0
        # GDS4522       0      NA       0       1       0       0       0       0       3       0       1
        # GDS4358       0       0      NA       5       0       0       0     263     410       0     599
        # GDS4218       0       0       0      NA       0       0       0      23      48       0     139
        # GDS4136       0       0       0       0      NA       0       0       0       0       0       0
        # GDS2821       0       0       0       1       0      NA       0       3       9       0      14
        # GDS1917       0       0       0       1       0       0      NA       0       0       0       0
        # GDS5204       0       0       0     479      19       1       5      NA    1578       0    2037
        # GDS4135       0       0       0       0       0       0       0       0      NA       0   10127
        # GDS2795       0       0       0       0       0       0       0       3       0      NA       0
        # GDS1962       1       0       0     997      21       3       3    2989       0       1      NA
        
        library(ggnetwork)
        out <- NULL
        for( src in rownames(deg_matrix) ) {
            for( tgt in colnames(deg_matrix) ) {
                arow <- data_frame( src = src, tgt = tgt, ndeg = deg_matrix[src, tgt] )
                if( is.null(out) ) {
                    out = arow
                } else {
                    out = bind_rows(out, arow)
                }
            }
        }
        
        library(network)
        library(sna)
        
        my_deg <- sapply(all_sam, function(sam) sam$siggenes.table$ngenes.up + sam$siggenes.table$ngenes.lo )
    
        dnet <- network(deg_matrix, #vertex.attr = my_deg, vertex.attrnames = rep('thisDEG',length(my_deg)),
                        directed = FALSE, bipartite = FALSE)
        dnet <- set.edge.value(dnet, 'n_DEG_lo',   deg_matrix)
        dnet <- set.edge.value(dnet, 'n_DEG_up', t(deg_matrix))
        as.matrix(dnet,attrname='n_DEG_lo')
        as.matrix(dnet,attrname='n_DEG_up')
        
        
        library(ggnetwork)
        
        set.seed(1234)
        png("caprice/RESULT/sam-deg-network.png", width = 960, height=960)
            ggplot(dnet, aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_edges(color = "grey50", aes(size=n_DEG_up+n_DEG_lo) ) +
            geom_edgetext_repel(aes(label = n_DEG_up), color = "black", size = 7,
                                fill = "#E6AAAA", box.padding = unit(1.5, "lines") ) +
            geom_edgetext_repel(aes(label = n_DEG_lo), color = "black", size = 7,
                                fill = "skyblue", box.padding = unit(.5, "lines")) +
            geom_nodes(color = "white", size = 17) +
            geom_nodetext(aes(label = paste0("GDS\n",substr(vertex.names,4,100)) ),
                          fontface = "bold", color="black", size=8, lineheight=.8) +
            theme_blank()
        dev.off()
        
    } else {
        write.table(deg_up[[i]], paste(datasets_num[[i]], '_genes_up.csv', sep=''), sep=',')
        write.table(deg_lo[[i]], paste(datasets_num[[i]], '_genes_lo.csv', sep=''), sep=',')         
    }
}


# Subset expr_vals to just DEG genes
exp_deg = datasets[[1]][which(rownames(datasets[[i]]) %in% deg_lo[[1]]),]


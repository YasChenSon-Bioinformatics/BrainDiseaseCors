source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library('samr')
library('limma')

source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

datasets_num = c(4523, 4522, 4358, 4218, 4136,
                 2821, 1917, 5204, 4135, 2795,
                 1962)

GDSl <- download_GDSs(GDSnumberv = datasets_num,
                      skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
datasets <- extractMatrixFromEset(ESETl)

# COMPUTE SAM
#
# SUPER HEAVY ( ~ 17 GB). Use a generated csv file below instead
#
bool_to_num = function(x){ if (x==TRUE) 1 else 2}
all_sam = list()

for (i in seq_along(datasets)) {
    df <- datasets[[i]]
    expCond <- GDSl[[i]]@dataTable@columns
    
    message(datasets_num[[i]], ' -----------------------------------------------')
    
    # if (i < 8) { next; } # uncomment to fast forward
    #if (i == 4) {

    # Possible levels used for control are:
    # 
    # "control"
    # "normal" (2795)
    # "young" (aging study 5204)
    # "Braak stage I-II" (alz 4135)
    # "non-tumor" (scf/angiogenesis 1962)
    
    control_strings = c("control", "normal", "Braak_stage_I-II", "non-tumor")
    dz_ctrl_boolv = grepl(paste(control_strings, collapse = "|"), expCond$disease.state)
    
    # Aging study has no disease.state column - use "age"
    # Note - grepl doesn't seem to like special characters
    if (datasets_num[[i]] == 5204) {
        age_control_strings = c("young", "middle")
        dz_ctrl_boolv = grepl(paste(age_control_strings, collapse = "|"), expCond$age)
    }
    
    y = apply(as.data.frame(dz_ctrl_boolv), 1, bool_to_num)
    
    # FIXME: change nperm from 2 (debug) to 1000 (production)
    samfit <- SAM(df, y, resp.type="Two class unpaired", nperms=100, fdr.output = 0.1,
                  geneid = rownames(df) )
    all_sam[[i]] = samfit
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
        # GDS4522       0      NA       0       0       0       0       0       0       0       0       0
        # GDS4358       0       0      NA       0       0       0       0       0       0       0       0
        # GDS4218       0       1       0      NA       0       0       0       0       0       0       0
        # GDS4136       0       0       0       0      NA       0       0       0       0       0       0
        # GDS2821       0       0       0       1       0      NA       0       0       0       0       0
        # GDS1917       0       0       0       1       0       0      NA       0       0       0       0
        # GDS5204       0       1       0     564       0       2       2      NA       0       0       0
        # GDS4135       0       0       0       0       0       0       0       0      NA       0       0
        # GDS2795       0       0       0       0       0       0       0       0       0      NA       0
        # GDS1962       0       1       0    1054       0       5       2    4007       0       0      NA
        
        library(ggnetwork)

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


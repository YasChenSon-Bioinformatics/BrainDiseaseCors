
source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')

setGlobalConstantList()
loadLibraries()

GDSlist <- download_GDSs(GDSnumberv = c( 5204, 4838, 4532, 4523, 4522,
                                         4477, 4358, 4231, 4218, 4136,
                                         4135, 3502, 2821, 2795, 2154,
                                         1962, 1917),
                         skipv = c('not-GPL570'),
                         suppress = FALSE)
ESETlist <- convertGDS2ESET( GDSlist )
Mlist <- extractMatrixFromEset( ESETlist )

sapply(seq_along(ESETlist), function(i) {table(is.na(exprs(ESETlist[[i]])) + 0) })
       
library(pheatmap)

binary_matrix <- binary_matrix[1:1000, 1:5]

i <- 5 # for test purpose
for ( i in seq_along(ESETlist) ) {
    thisESET <- ESETlist[[i]]
    GDS_number <- thisESET@experimentData@other$dataset_id[1]
    binary_matrix <- t(is.na(exprs(thisESET)) + 0) # transpose for wider macbook screen width
    sorted_matrix <- binary_matrix[ names(sort(rowSums(binary_matrix))),  ]
    
    n_NA <- sum(binary_matrix)
    n_sample <- nrow(sorted_matrix)
    png(paste0("EDA/NA/heatmap-", GDS_number, ".png"), width=2560, height=1600 ) # pixel
    if( n_NA > 0 ) {
        pheatmap(sorted_matrix,
                 cluster_row = FALSE, # faster
                 cluster_col = FALSE, # faster
                 color=gray.colors(2,start=1,end=0), # black: NA. white : otherwise
                 show_colnames = FALSE, # too cluttered
                 main = paste0("NA locations of ", GDS_number,
                               " (Number of NAs : ", n_NA, " ) / (Number of Samples : ", n_sample, " )",
                               " =  about ", round(n_NA / n_sample,0) , " probes are missing" ))
    } else { # no NA in datasets
        plot(0, type='n', main = paste0("NA locations of ", GDS_number))
        text(1,0,'No NA in datasets. Instead show the original Expression Values')
        pheatmap(t(exprs(thisESET)),
                 cluster_row = FALSE, # faster
                 cluster_col = FALSE, # faster
                 show_colnames = FALSE, # too cluttered
                 main = paste0("Expression Values of ", GDS_number))
    }
    dev.off()
    message('Plotted ', i)
}

# Then concatenated by 'convert *.png 17datasets-heatmap.pdf' on terminal
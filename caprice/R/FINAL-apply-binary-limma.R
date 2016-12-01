source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
Ml <- extractMatrixFromEset(ESETl)
topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=400, p_threshold=.1, method='fdr')
#topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, p_threshold = .1, method = 'fdr')
deg_matrix <- build_deg_matrix(topped)

deg_matrix # nTopGene= 400
#         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# GDS5204      NA       0       0       0       0       1       0       0       1       1       0
# GDS4522       0      NA      11       4       1       0       1       5       1       3       0
# GDS4523       1       5      NA       1       2       0       2       5       2       0       0
# GDS4358       0       0       0      NA       0       2       2       2       3       2       0
# GDS4218       0       5       0       0      NA       0       0       0       1       1       0
# GDS4136       5       0       0       0       0      NA       5       1       4       3       3
# GDS4135       0       0       1       0       0       0      NA       1       5       3       5
# GDS2821       4       2       1       0       0       1       0      NA       1       2       1
# GDS2795       4       0       0       0       0       0       0       1      NA       2       1
# GDS1962      13       0       2       0       0       3       0       0       1      NA       1
# GDS1917       0       1       1       0       0       1       0       0       0       0      NA

library(extrafont)
loadfonts()

library(ggnetwork)
library(network)
library(sna)

fixed_deg_layout <-
    matrix(data = c(  1,   .45,  # 5204
                      .05, .5 ,  # 4522
                      .2 , .45,  # 4523
                      0, .9,     # 4358
                      .2,  0,    # 4218
                      .75,  1,  # 4136
                      .4 , .85 ,  # 4135
                     .85,   .63,  # 2821
                      .6,   .35, # 2795
                      .4,   .68, # 1962
                      .1,   1    # 1917
    ),
    ncol=2, byrow = TRUE)

dzsummary <-
    list(
        GDS5204 = '5204\nAging',
        GDS4522 = '4522\nSchizophrenia', # 'Schizophrenia',
        GDS4523 = '4523\nSchizophrenia',
        GDS4358 = '4358\nHIV',
        GDS4218 = '4218\nMultiple\nSclerosis',
        GDS4136 = '4136\nAlzheimer',#'er',
        GDS4135 = '4135\nBraak',
        GDS2821 = '2821\nParkinson',#, 'son',
        GDS2795 = '2795\nNeurofibrillary\nTangle', #'fibrillary_Tangle',
        GDS1962 = '1962\nTumor',
        GDS1917 = '1917\nSchizophrenia'
    )

pdf("caprice/RESULT/binary-limma-deg-network.pdf", width = 12, height=12)
for( nn in seq(100,500, by=25) ){
    topped <- applyTtestToGeneExpressionMatrices(
                    Ml, GDSl, nTopGene=nn, p_threshold=.1, method='fdr'
                )
    deg_matrix <- build_deg_matrix(topped)
    
    dnet <- network(deg_matrix,
                    directed = FALSE, bipartite = FALSE)
    dnet <- set.edge.value(dnet, 'n_DEG_lo',   deg_matrix)
    dnet <- set.edge.value(dnet, 'n_DEG_up', t(deg_matrix))
    dnet <- set.edge.value(dnet, 'n_DEG',
                           matrix(paste0(t(deg_matrix), '\n', deg_matrix),
                                  nrow=nrow(deg_matrix)))
    dnet %v% "gds_dz" <- unlist(dzsummary[ rownames(deg_matrix) ])
    dnet %v% "dz" <- unlist(dzsummary[ rownames(deg_matrix) ]) %>% gsub('[0-9]+\n','',.)
    
    set.seed(1234)
    p <-
        ggplot(ggnetwork(dnet, layout = fixed_deg_layout),
               aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size=n_DEG_up+n_DEG_lo,
                       alpha=ifelse(n_DEG_up+n_DEG_lo<2,0,
                                    ifelse(n_DEG_up+n_DEG_lo>5,1,.5))),
                   show.legend = FALSE) +
        geom_edgetext(aes(label = n_DEG), color = "black", size = 7,
                      fill = "#E2E9EC", lineheight=.8) +
        geom_nodes(aes(color = dz), size = 27, shape=15) +
        #geom_nodes(color = "white", size = 24) +
        geom_nodetext(aes(label = gds_dz),
                      fontface = "bold", color="white", size=7, lineheight=.8) +
        theme_blank() + labs(title = paste0("Number of DEGs for each dataset = ", nn))
    print(p)
}
dev.off()

final_out <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=400,
                                                p_threshold=.05, method='fdr', type ='p')

sapply(final_out, function(x) x$table)
sapply(final_out, function(x) ifelse(nrow(x$table) > 0, x$gds, "") )

final_deg_matrix <- build_deg_matrix(final_out)
# all zero matrix. this is understandable considering the correction method is
# Benjamini and Hochberg (1995), which is mentioned in the SAM paper (2001).

final_out <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=300,
                                                p_threshold=.05, method='fdr')
final_deg_matrix <- build_deg_matrix(final_out)

png("caprice/RESULT/binary-limma-deg-network.png", width = 960, height=960)
    dnet <- network(final_deg_matrix,
                    directed = FALSE, bipartite = FALSE)
    dnet <- set.edge.value(dnet, 'n_DEG_lo',   final_deg_matrix)
    dnet <- set.edge.value(dnet, 'n_DEG_up', t(final_deg_matrix))
    dnet <- set.edge.value(dnet, 'n_DEG',
                           matrix(paste0(t(final_deg_matrix), '\n', final_deg_matrix),
                                  nrow=nrow(final_deg_matrix)))
    dnet %v% "gds_dz" <- unlist(dzsummary[ rownames(final_deg_matrix) ]) %>% substr(., 1,10)
    dnet %v% "disease" <- unlist(dzsummary[ rownames(final_deg_matrix) ]) %>% gsub('[0-9]+\n','',.)
    
    set.seed(1234)
    p <-
        ggplot(ggnetwork(dnet, layout = fixed_deg_layout),
               aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size=n_DEG_up+n_DEG_lo,
                       alpha=ifelse(n_DEG_up+n_DEG_lo<2,0,
                             ifelse(n_DEG_up+n_DEG_lo>5,1,.5))),
                   show.legend = FALSE) +
        geom_edgetext(aes(label = n_DEG), color = "black", size = 7,
                      fill = "#E2E9EC", lineheight=.8) +
        geom_nodes(aes(color = disease), size = 27, shape=15) +
        #geom_nodes(color = "white", size = 24) +
        geom_nodetext(aes(label = gds_dz),
                      fontface = "bold", color="white", size=7, lineheight=.8) +
        theme_blank() + labs(title = paste0("Number of DEGs for each dataset = ", nn),
                             subtitle = "First/Second Number in Edges denote the Number of Up-/Down- regurated probes") + theme(legend.text=element_text(size=18, face="bold"), legend.key.size = unit(.1, "cm"))
    print(p)
dev.off()

write.csv(attr(final_deg_matrix, 'df'),
          file='caprice/RESULT/11GDS_binary_limma_results.csv', quote=FALSE, row.names=FALSE)

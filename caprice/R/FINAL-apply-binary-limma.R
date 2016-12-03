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

# Affymetrix probe id -> uniplot
library(hgu133plus2.db)


# Source: http://www.reactome.org/download/current/UniProt2Reactome.txt
#    and: http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
#
# The above URL's UniProt2Reactome.txt is more than 30 MB.
#                (UniProt2Reactome_All_Levels.txt is more than 130 MB)
# Thus I extracted necessary parts ( < 4MB and 10MB each ) by the following commands 
#
# read.delim("http://www.reactome.org/download/current/UniProt2Reactome.txt", sep='\t',
#            header = FALSE, col.names = c('uni', 'id', 'url', 'pathway', 'type', 'sp'), 
#            stringsAsFactors = FALSE)  %>%
# filter(sp == 'Homo sapiens') %>% dplyr::select(uni, pathway, url) %>%
# group_by(uni, pathway) %>% dplyr::slice(1) %>% group_by() %>%
# write.table(., file = "caprice/MAP/UniProt2Reactome.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#
# read.delim("http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt", sep='\t',
#            header = FALSE, col.names = c('uni', 'id', 'url', 'pathway', 'type', 'sp'),
#            stringsAsFactors = FALSE)  %>%
#    filter(sp == 'Homo sapiens') %>% dplyr::select(uni, pathway, url) %>%
#    group_by(uni, pathway) %>% dplyr::slice(1) %>% group_by() %>%
#    write.table(., file = "caprice/MAP/UniProt2Reactome_All_Levels.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#
# group_by(uni, pathway) %>% dplyr::slice(1) is because uni-pathway is duplicated for some genes for some reason 
#
# Note: There are SEVERAL pathway databases. We use REACTOME database here, which we learend in classes.
# See: http://www.mpb.unige.ch/reports/rap_Jia_Li.pdf
#

gene2pathwaydf <-
    #read.delim("caprice/MAP/UniProt2Reactome_All_Levels.tsv",  # too minute
    read.delim("caprice/MAP/UniProt2Reactome.tsv",
               sep='\t', header = TRUE, stringsAsFactors=FALSE) %>%
    mutate( pathway = gsub(' ','_',pathway) )

# do pathway enrichment analysis

# final_out (nTopGene=300) returns too few common pathways.
topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=500,
                                             p_threshold=.05, method='fdr')

peal <- list() # pea list
for( i in seq_along(gdsv) ){
    message("  ", gdsv[i])
    peal[[i]] <- do_PathwayEnrichmentAnalysis(gene2pathwaydf, topped, gdsv, gds = gdsv[i])
}

pea_matrix <- do.call(cbind, lapply( peal, function(x) x$pval ))
colnames(pea_matrix) <- gdsv
rownames(pea_matrix) <- peal[[1]]$pathway

as.data.frame(pea_matrix) %>% rownames_to_column("pathway") %>%
    dplyr::select(pathway, GDS5204, GDS1962) %>% filter( GDS5204 < .01 & GDS1962 < .05 )

#
# The following pathways should be related to both diseases.
#
#                          pathway                        GDS5204_pval GDS1962_pval
# 1                                          Ca2+_pathway 1.607336e-03 4.039154e-02
# 2                                  Cam-PDE_1_activation 2.682637e-03 2.801654e-03
# 3  CREB_phosphorylation_through_the_activation_of_CaMKK 2.682637e-03 2.801654e-03
# 4 Interactions_of_neurexins_and_neuroligins_at_synapses 1.754454e-03 4.599099e-05
# 5                 MASTL_Facilitates_Mitotic_Progression 1.051869e-03 1.926011e-02
# 6                        Phase_0_-_rapid_depolarisation 3.729520e-04 1.673485e-02
# 7                         Trafficking_of_AMPA_receptors 2.199996e-05 5.670350e-03

arbitrary_p_threshold <- .01 # Just for visualization

pea_bitmatrix <- pea_matrix < arbitrary_p_threshold 

#  [ 11 x 1740 ]       [ 1740 x 11 ]
common_pathway_matrix <- t(pea_bitmatrix) %*% pea_bitmatrix
#     diagonal elements of common_pathway_matrix are the numbers of < arbitrary_p_threshold.
#              can be confirmed as: colSums(pea_bitmatrix) 
# non-diagonal elements of common_pathway_matrix are the numbers of common pathways between two GDSs.
common_pathway_matrix
#         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# GDS5204      52       0       0       0       4       4       0       1       0      11       0
# GDS4522       0       3       0       1       0       0       0       0       0       0       1
# GDS4523       0       0      10       2       0       0       0       0       0       0       0
# GDS4358       0       1       2      51       0       0       1       2       5       0       1
# GDS4218       4       0       0       0      13       0       0       0       0       1       0
# GDS4136       4       0       0       0       0       9       0       0       0       2       0
# GDS4135       0       0       0       1       0       0       6       0       1       0       0
# GDS2821       1       0       0       2       0       0       0      11       0       1       0
# GDS2795       0       0       0       5       0       0       1       0      31      12       0
# GDS1962      11       0       0       0       1       2       0       1      12      73       0
# GDS1917       0       1       0       1       0       0       0       0       0       0      18


cnet <- network(common_pathway_matrix,
                directed = FALSE, bipartite = FALSE)
cnet <- set.edge.value(cnet, 'n_common',   common_pathway_matrix)
cnet %v% "gds_dz" <- unlist(dzsummary[ rownames(common_pathway_matrix) ]) %>% substr(.,1,10)
cnet %v% "disease" <- unlist(dzsummary[ rownames(common_pathway_matrix) ]) %>% gsub('[0-9]+\n','',.)


ggplot(ggnetwork(cnet, layout = fixed_deg_layout),
       aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size=n_common,
                   alpha=ifelse(n_common<2,0,
                                ifelse(n_common>10,1,.5))),
               show.legend = FALSE) +
    geom_edgetext(aes(label = n_common), color = "black", size = 7,
                  fill = "#E2E9EC", lineheight=.8) +
    geom_nodes(aes(color = disease), size = 27, shape=15) +
    #geom_nodes(color = "white", size = 24) +
    geom_nodetext(aes(label = gds_dz),
                  fontface = "bold", color="white", size=7, lineheight=.8) +
    theme_blank() + labs(title = paste0("Number of DEGs for each dataset = ", 500),
                         subtitle = "Numbers in edges denote the number of common pathways") + theme(legend.text=element_text(size=18, face="bold"), legend.key.size = unit(.1, "cm"))



  
 
# uni                                               pathway                                                url n_path n_enriched.x       pval.x n_enriched.y     pval.y
# <chr>                                                 <chr>                                              <chr>  <dbl>        <dbl>        <dbl>        <dbl>      <dbl>
# P62158                                 Activation_of_CaMK_IV  http://reactome.org/PathwayBrowser/#/R-HSA-442745      4            4 2.065123e-05            2 0.02714004
# P54750                                  Cam-PDE_1_activation  http://reactome.org/PathwayBrowser/#/R-HSA-111957      4            2 2.494469e-02            2 0.02714004
# P16220              CaMK_IV-mediated_phosphorylation_of_CREB  http://reactome.org/PathwayBrowser/#/R-HSA-111932      5            4 9.770601e-05            2 0.04313568
# O43865             CLEC7A_(Dectin-1)_induces_NFAT_activation http://reactome.org/PathwayBrowser/#/R-HSA-5607763     12            3 4.266570e-02            3 0.04771767
# P16220  CREB_phosphorylation_through_the_activation_of_CaMKK  http://reactome.org/PathwayBrowser/#/R-HSA-442717      4            2 2.494469e-02            2 0.02714004
# O14490 Interactions_of_neurexins_and_neuroligins_at_synapses http://reactome.org/PathwayBrowser/#/R-HSA-6794361     87           11 3.179326e-02           10 0.08547570
# O43768                 MASTL_Facilitates_Mitotic_Progression http://reactome.org/PathwayBrowser/#/R-HSA-2465910     13            5 1.130047e-03            3 0.05887312
    
# Useful:
# CAMK IV in epithelial ovarian cancer. https://www.ncbi.nlm.nih.gov/pubmed/12065094
# https://www.ncbi.nlm.nih.gov/pubmed/8194751
# https://www.ncbi.nlm.nih.gov/pubmed/15680915


    
# AnnotationDbi::select(hgu133plus2.db, keys="241672_at", columns="UNIPROT")
# gene2pathwaydf %>% filter( uni == 'A2A2V5' ) # secondary Q8N469

# nDEG = 1000 for each dataset does not show so different results.
# Checked by:
# 
# n1000 <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=1000,
#                                                 p_threshold=.05, method='fdr')
# n1000_deg_matrix <- build_deg_matrix(n1000)
# 
# TBL5204 <- n1000[[ which(gdsv == 'GDS5204') ]]$table
# TBL1962 <- n1000[[ which(gdsv == 'GDS1962') ]]$table
# 
# PEA5204 <- do_pea(gene2pathwaydf, rownames(TBL5204), p_threshold = .2)
# PEA1962 <- do_pea(gene2pathwaydf, rownames(TBL1962), p_threshold = .2, nopath = TRUE)
# 
# joined <- inner_join(PEA5204, PEA1962, by='pathway')
# gene2pathwaydf %>% filter( pathway %in% joined$pathway ) %>% group_by(pathway) %>% dplyr::slice(1) %>%
#     group_by() %>% left_join(., joined, by="pathway")
# 




#AnnotationDbi::select(hgu133plus2.db, keys="241672_at", columns="ENSEMBL")
# http://www.genome.jp/dbget-bin/www_bget?hsa:400120
# no pathway reported ???

for( gds in gdsv ){
    print(do_pea(rownames(topped[[ which(gdsv == gds) ]]$table)))
}
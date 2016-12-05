source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
Ml <- extractMatrixFromEset(ESETl)

# Analysis 1: Differentially Expressed Genes (DEG) by Linear Model (binary)

topped <- applyTtestToGeneExpressionMatrices(Ml,GDSl,nTopGene=500,method='fdr',removeGenderEffect=TRUE)
#topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, p_threshold = .1, method = 'fdr')
deg_matrix <- build_deg_matrix(topped)

deg_matrix # nTopGene= 500
# Upper-triangles: up-regulated
# Lower-triangles: down-regulated
#
#         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# GDS5204      NA       0       0       2       0       1       0       0       1       1       0
# GDS4522       1      NA      15       4       1       1       2       7       1       3       1
# GDS4523       1       7      NA       4       3       1       4       6       3       2       0
# GDS4358       0       0       0      NA       0       6       4       0       5       4       2
# GDS4218       0      10       0       0      NA       0       0       1       1       1       1
# GDS4136       9       1       0       0       0      NA       7       2       4       3       3
# GDS4135       0       0       1       0       0       0      NA       3       9       3       7
# GDS2821       6       2       2       0       2       2       0      NA       1       3       1
# GDS2795       5       0       0       0       0       0       0       1      NA       3       1
# GDS1962      22       0       2       0       0       3       1       0       1      NA       2
# GDS1917       1       1       1       0       1       1       0       1       1       0      NA
write.csv(deg_matrix, file = "caprice/RESULT/binary-limma-1-deg-matrix.csv")

# Affymetrix probe id -> uniplot
library(hgu133plus2.db)
probe2unidf <- 
AnnotationDbi::select(hgu133plus2.db,
                      keys=intersect(rownames(topped[[1]]$table), rownames(topped[[10]]$table)), # 23 probes
                      columns="UNIPROT") %>% group_by(PROBEID) %>%
    summarize( uni = paste0(UNIPROT, collapse='|') ) 

gds5204info <- GDSl[[1]]@dataTable@columns %>% dplyr::select(sample, age) %>%
    mutate( disease = buildDesignMatrix(GDSl[[1]]) %>% as.data.frame %>% .$dz_TRUE %>% as.factor) %>%
    dplyr::rename( original = age ) %>% mutate( gds = 'GDS5204' )
gds1962info <- GDSl[[10]]@dataTable@columns %>% dplyr::select(sample, disease.state) %>%
    mutate( disease = buildDesignMatrix(GDSl[[10]]) %>% as.data.frame %>% .$dz_TRUE %>% as.factor) %>%
    dplyr::rename( original = disease.state )  %>% mutate( gds = 'GDS1962' )

gene_expression_values <-
    cbind(exprs(ESETl[[1]]), exprs(ESETl[[10]]))[ probe2unidf$PROBEID, ] %>%
    as.data.frame %>% # dplyr cannot handle a class 'matrix'
    rownames_to_column('probe') %>%
    gather(key=smpl, value=eval, -probe) %>%
    left_join(. , bind_rows(gds5204info, gds1962info), c('smpl' = 'sample')) %>%
    group_by(probe, disease, gds) %>% mutate( hypothesis = mean(eval) ) %>%
    group_by(gds) %>% mutate( r = row_number() ) %>%
    group_by(gds, disease) %>% mutate( x = min(r), xend = max(r) ) %>%
    left_join(., probe2unidf, by=c('probe' = 'PROBEID')) %>%
    left_join(
        .,
        bind_rows(topped[[1]]$table %>% rownames_to_column('probe') %>% mutate( gds = 'GDS5204' ),
                  topped[[10]]$table %>% rownames_to_column('probe') %>% mutate( gds = 'GDS1962' )),
        by = c('gds', 'probe')
    ) %>%
    mutate( lab = paste0(probe, ' (', uni, ')\nt = ', round(t,4), '\n(adjusted-p = ', round(adj.P.Val,10), ')' ) )
    
mycap <- 'source code: FINAL-apply-binary-limma.R in github.com/YasChenSon-Bioinformatics/BrainDiseaseCors'

png('caprice/RESULT/binary-limma-2-23probes-GDS5204.png', height=960, width=1200)
    gene_expression_values %>% filter( gds == 'GDS5204' ) %>%
    ggplot() +
    geom_point(aes(x=r, y=eval, color=disease)) +
    geom_segment(aes(x=x, xend=xend, y=hypothesis, yend=hypothesis, color=disease)) + # For Male
    facet_wrap( ~ lab, ncol =4 ) + # create subplots for each probe
    labs(title='Raw Expression Values of 23 Probes in GDS5204',
         subtitle = 'x-axis: Sample Index (hidden)   y-axis: Expression Value    Color: Disease or Not    Horizontal Bar: Hypothesis',
         caption = mycap,
         x = '', y = 'Expression Values') +
    theme( axis.text.x = element_blank(), axis.ticks.x = element_blank(),
           legend.position = c(.88, .1), text = element_text(size=14, face = 'bold') ) +
           scale_color_discrete(labels = c('Not Disease (<=70 years old)', 'Disease (>70 years old)')) +
    scale_y_continuous(breaks=seq(1,3.6,by=.2))
dev.off()

png('caprice/RESULT/binary-limma-2-23probes-GDS1962.png', height=960, width=1200)
gene_expression_values %>% filter( gds == 'GDS1962' ) %>%
    ggplot() +
    geom_point(aes(x=r, y=eval, color=disease), size = 1.2) +
    geom_segment(aes(x=x, xend=xend, y=hypothesis, yend=hypothesis, color=disease)) + # For Male
    facet_wrap( ~ lab, ncol = 4 ) + # create subplots for each probe
    labs(title='Raw Expression Values of 23 Probes in GDS1962',
         subtitle = 'x-axis: Sample Index (hidden)   y-axis: Expression Value    Color: Disease or Not    Horizontal Bar: Hypothesis',
         caption = mycap,
         x = '', y = 'Expression Values') +
    theme( axis.text.x = element_blank(), axis.ticks.x = element_blank(),
           legend.position = c(.88, .1), text = element_text(size=14, face = 'bold') ) +
    scale_color_discrete(labels = c('Non-Tumor', 'Tumor')) +
    scale_y_continuous(breaks=seq(2,15,by=1))
dev.off()

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
                      .35,   .68, # 1962
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

final_out <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=500,
                                                p_threshold=.05, method='fdr')
final_deg_matrix <- build_deg_matrix(final_out)

png("caprice/RESULT/binary-limma-3-deg-network.png", width = 960, height=960)
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
                       alpha=ifelse(n_DEG_up+n_DEG_lo<5,0,
                             ifelse(n_DEG_up+n_DEG_lo>15,1,.5))),
                   show.legend = FALSE) +
        geom_edgetext(aes(label = n_DEG
                          # size = ifelse(n_DEG_up+n_DEG_lo<5,1,
                          #               ifelse(n_DEG_up+n_DEG_lo>15,20,4))  # doesn't work
                              ), color = "black", size = 7,
                      fill = "#E2E9EC", lineheight=.8) +
        geom_nodes(aes(color = disease), size = 27, shape=15) +
        #geom_nodes(color = "white", size = 24) +
        geom_nodetext(aes(label = gds_dz),
                      fontface = "bold", color="white", size=7, lineheight=.8) +
        theme_blank() + labs(title = paste0("Commonly Differentially Expressed Probes / Number of DEGs for each dataset = ", nn),
                             subtitle = "First/Second Number in Edges denote the Number of Up-/Down- regurated probes",
                             caption = mycap) +
        theme(legend.text=element_text(size=18, face="bold"), legend.key.size = unit(.1, "cm"),
              text = element_text(face='bold'))
    print(p)
dev.off()

write.csv(attr(final_deg_matrix, 'df'),
          file='caprice/RESULT/11GDS_binary_limma_results.csv', quote=FALSE, row.names=FALSE)




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

# Note: (adjusted) p < .05 cutoff returns 10 probes for both GDS5204 and GDS1962 (no common probes).
# example:
#   topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=500,
#                                                  p_threshold=.05, method='fdr', type = 'p')
# This is too stringent and provides no information, so instead we return top 500 probes here.


peal <- list() # pea list
for( i in seq_along(gdsv) ){
    message("  ", gdsv[i])
    peal[[i]] <- do_PathwayEnrichmentAnalysis(gene2pathwaydf, topped, gdsv, gds = gdsv[i])
}

pea_matrix <- do.call(cbind, lapply( peal, function(x) x$pval ))
colnames(pea_matrix) <- gdsv
rownames(pea_matrix) <- peal[[1]]$pathway

as.data.frame(pea_matrix) %>% rownames_to_column("pathway") %>%
    dplyr::select(pathway, GDS5204, GDS1962) %>% filter( GDS5204 < .05 & GDS1962 < .05 )

#
# The following pathways should be related to both diseases.
#                                               pathway         GDS5204        GDS1962
# 1                                           Ca2+_pathway 1.607335614e-03 0.040391538103
# 2                                   Cam-PDE_1_activation 2.682637108e-03 0.002801654141
# 3              CLEC7A_(Dectin-1)_induces_NFAT_activation 2.226682269e-02 0.023203121181
# 4   CREB_phosphorylation_through_the_activation_of_CaMKK 2.682637108e-03 0.002801654141
# 5                                        DARPP-32_events 1.421792226e-02 0.015078244356
# 6                                     DSCAM_interactions 2.226682269e-02 0.023203121181
# 7  Interactions_of_neurexins_and_neuroligins_at_synapses 1.754453963e-03 0.000045990988
# 8                                  LGI-ADAM_interactions 3.532034703e-02 0.003186419020
# 9                  MASTL_Facilitates_Mitotic_Progression 1.051869167e-03 0.019260107293
# 10                        Phase_0_-_rapid_depolarisation 3.729520433e-04 0.016734851346
# 11               Synthesis_of_IP3_and_IP4_in_the_cytosol 1.770038672e-02 0.018758429177
# 12                         Trafficking_of_AMPA_receptors 2.199995844e-05 0.005670349630
#
# more details are in the following csv file

tmp <- left_join(peal[[1]], peal[[10]], by = 'pathway') %>%
    filter( pval.x < .05 & pval.y < .05 ) %>%
    mutate_each(funs(round(.,6)), starts_with('pval'), starts_with('fdr')) %>%
    dplyr::rename( n_pathway = n_path.x ) %>% dplyr::select(-n_path.y) %>%
    left_join(.,
              as.data.frame(pea_matrix) %>% rownames_to_column('pathway') %>% mutate_each(funs(round(.,4)), -pathway),
              by='pathway') %>%
    arrange(pval.x) %>%
    dplyr::select(pathway, n_pathway, starts_with('n_enriched'), starts_with('pval'),
                  starts_with('fdr'), starts_with('gene'), starts_with('GDS')) %>%
    dplyr::select(-GDS5204, -GDS1962)
colnames(tmp) <- gsub('.x', '_GDS5204',colnames(tmp)) %>% gsub('\\.y', '_GDS1962', .) %>% gsub('n_enriched', 'n_DEG', .)

tmp <- gene2pathwaydf %>% filter( pathway %in% tmp$pathway ) %>%
    group_by(pathway) %>% dplyr::slice(1) %>%
    group_by() %>% dplyr::select(-uni) %>% left_join(tmp, ., by="pathway") %>%
    dplyr::rename( nGeneInPathway = n_pathway )
write.csv(tmp, file='caprice/RESULT/binary-limma-5-common-pathways-GDS5204-GDS1962.csv', row.names=FALSE, quote=FALSE)

arbitrary_p_threshold <- .05 # Just for visualization

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

png('caprice/RESULT/binary-limma-4-pathway-network.png', width=960, height=960)
    ggplot(ggnetwork(cnet, layout = fixed_deg_layout),
       aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size=n_common,
                   alpha=ifelse(n_common<2,0,
                                ifelse(n_common>10,1,.5))),
               show.legend = FALSE) +
    geom_edgelabel(aes(label = n_common), color = "black", size = 7,
                  fill = "#FFFFFF", lineheight=.8) +
    geom_nodes(aes(color = disease), size = 27, shape=15) +
    #geom_nodes(color = "white", size = 24) +
    geom_nodetext(aes(label = gds_dz),
                  fontface = "bold", color="white", size=7, lineheight=.8) +
    theme_blank() + labs(title = paste0("COMMON PATHWAYS / Number of DEGs for each dataset = ", 500),
                         subtitle = "Numbers in edges denote the number of common pathways",
                         caption = mycap) +
        theme(legend.text=element_text(size=18, face="bold"), legend.key.size = unit(.1, "cm"),
              text=element_text(face='bold'))
dev.off()
    
# Useful:
# CAMK IV in epithelial ovarian cancer. https://www.ncbi.nlm.nih.gov/pubmed/12065094
# https://www.ncbi.nlm.nih.gov/pubmed/8194751
# https://www.ncbi.nlm.nih.gov/pubmed/15680915


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

age_values <-
    exprs(ESETl[[1]])[ probe2unidf$PROBEID, ] %>%
    as.data.frame %>% # dplyr cannot handle a class 'matrix'
    rownames_to_column('probe') %>%
    gather(key=smpl, value=eval, -probe) %>%
    left_join(. , tmp %>% mutate( a = gsub(".* ([0-9\\.]+)\\s+years.*","\\1", description, perl = TRUE) ), c('smpl' = 'sample'))
    
age_values %>% ggplot() + geom_point(aes(x=age,y=eval)) + facet_wrap( ~ probe )
age_values %>% ggplot(aes(x=as.numeric(a),y=eval)) + geom_point() + geom_smooth(method="lm") + geom_vline(xintercept = 70, color="red")+ facet_wrap( ~ probe )
    
source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
Ml <- extractMatrixFromEset(ESETl)

GDS4522probev <- rownames(Ml[[ which(sapply(GDSl, function(x) x@header$dataset_id[1]) == 'GDS4522') ]])

Ml2 <- sapply(Ml, function(m) m[ rownames(m) %in% GDS4522probev, ] )

sapply(Ml2, function(x) dim(x))

topped_subset <- applyTtestToGeneExpressionMatrices(Ml2, GDSl, nTopGene = 500)

deg_subset_matrix <- build_deg_matrix(topped_subset)

deg_subset_matrix

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

dnet <- network(deg_subset_matrix,
                directed = FALSE, bipartite = FALSE)
dnet <- set.edge.value(dnet, 'n_DEG_lo',   deg_subset_matrix)
dnet <- set.edge.value(dnet, 'n_DEG_up', t(deg_subset_matrix))
dnet <- set.edge.value(dnet, 'n_DEG',
                       matrix(paste0(t(deg_subset_matrix), '\n', deg_subset_matrix),
                              nrow=nrow(deg_subset_matrix)))
dnet %v% "gds_dz" <- unlist(dzsummary[ rownames(deg_subset_matrix) ]) %>% substr(., 1,10)
dnet %v% "disease" <- unlist(dzsummary[ rownames(deg_subset_matrix) ]) %>% gsub('[0-9]+\n','',.)

set.seed(1234)
p <-
    ggplot(ggnetwork(dnet, layout = fixed_deg_layout),
           aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size=n_DEG_up+n_DEG_lo,
                   alpha=ifelse(n_DEG_up+n_DEG_lo<5,0,
                                ifelse(n_DEG_up+n_DEG_lo>25,1,.5))),
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
    theme_blank() + labs(title = paste0("Commonly Differentially Expressed Probes / Number of DEGs for each dataset = 500"),
                         subtitle = "First/Second Number in Edges denote the Number of Up-/Down- regurated probes"
                         ) +
    theme(legend.text=element_text(size=18, face="bold"), legend.key.size = unit(.1, "cm"),
          text = element_text(face='bold'))
print(p)
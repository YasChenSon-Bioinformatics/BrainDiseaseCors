degdf <- read.csv('caprice/RESULT/11GDS_sam_results.csv', stringsAsFactors = FALSE) %>%
    mutate( gds = paste0('GDS', gds) )

## Deg matrix
# FIXME: duplicate build_deg_marix()

type_up <- 'up'
type_lo <- 'lo'

gdsnov <- sapply(GDSl, function(x) x@header$dataset_id[1])  # assume length(GDSl) is 11

k <- 1
i <- 2

deg_matrix <- matrix( NA, ncol = length(gdsnov), nrow = length(gdsnov) )

rownames(deg_matrix) <- gdsnov
colnames(deg_matrix) <- gdsnov

for( i in seq_along(gdsnov) ){
    for( k in i:length(gdsnov) ){
        # AveExpr == 0 is a placeholder used in applyTtestToGeneExpressionMatrices
        if( nrow( degdf %>% filter( gds == gdsnov[k] ) ) == 0 ||
            nrow( degdf %>% filter( gds == gdsnov[i] ) ) == 0 ) {
            deg_matrix[i,k] <- 0
            deg_matrix[k,i] <- 0
            next
        }
        
        gds_i <- gdsnov[i]
        gds_k <- gdsnov[k]
        
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




### pathway enrichment analysis

gene2pathwaydf <-
    read.delim("caprice/MAP/UniProt2Reactome.tsv",
               sep='\t', header = TRUE, stringsAsFactors=FALSE) %>%
    mutate( pathway = gsub(' ','_',pathway) )


sampeal <- list() # pea list
peaed_gds <- c()
for( i in seq_along(gdsv) ){ # FIXME gdsv
    message("  ", gdsv[i])
    probev <- univ %>% filter( gds == gdsv[i] ) %>% .$probe
    if (length(probev) < 3){
        message("too few probes") # FIXME 
        next
    } else {
        peaed_gds <- c(peaed_gds, gdsv[i])
    }
    sampeal[[i]] <- do_pea(gene2pathwaydf, probev, p_threshold = 1, nopath=TRUE)
}

pea_matrix <- do.call(cbind, lapply( sampeal, function(x) x$pval ))
colnames(pea_matrix) <- peaed_gds
rownames(pea_matrix) <- sampeal[[1]]$pathway

as.data.frame(pea_matrix) %>% rownames_to_column("pathway") %>%
    dplyr::select(pathway, GDS5204, GDS1962) %>% filter( GDS5204 < .05 & GDS1962 < .05 )

# pathway     GDS5204    GDS1962
# 1  Activation_of_the_AP-1_family_of_transcription_factors 0.009516941 0.03109728
# 2    Advanced_glycosylation_endproduct_receptor_signaling 0.012356084 0.04024217
# 3      CREB_phosphorylation_through_the_activation_of_Ras 0.012356084 0.04024217
# 4                                        CYP2E1_reactions 0.010464137 0.03415493
# 5                                    ERKs_are_inactivated 0.011410518 0.03720322
# 6        Gastrin-CREB_signalling_pathway_via_PKC_and_MAPK 0.008568929 0.02803025
# 7     Golgi_Cisternae_Pericentriolar_Stack_Reorganization 0.013300836 0.04327183
# 8                                 MAPK1_(ERK2)_activation 0.009516941 0.03109728
# 9            Negative_feedback_regulation_of_MAPK_pathway 0.005719991 0.01877257
# 10                     RAF-independent_MAPK1/3_activation 0.012356084 0.04024217
# 11                                         RSK_activation 0.005719991 0.01877257
# 12                                     Signal_attenuation 0.009516941 0.03109728
# 13                       Spry_regulation_of_FGF_signaling 0.015187901 0.04930334
# 14        Zinc_influx_into_cells_by_the_SLC39_gene_family 0.009516941 0.03109728


as.data.frame(pea_matrix) %>% rownames_to_column("pathway") %>%
     dplyr::select(pathway, GDS5204, GDS1962) %>% filter( grepl("CaMK|Ca2+",pathway) ) %>% arrange(GDS5204)
#                                                               pathway GDS5204 GDS1962
# 1                                          Ca2+_activated_K+_channels       1       1
# 2                                                        Ca2+_pathway       1       1
# 3                            CaMK_IV-mediated_phosphorylation_of_CREB       1       1
# 4               CREB_phosphorylation_through_the_activation_of_CaMKII       1       1
# 5                CREB_phosphorylation_through_the_activation_of_CaMKK       1       1
# 6                                  Elevation_of_cytosolic_Ca2+_levels       1       1
# 7 Inhibition__of_voltage_gated_Ca2+_channels_via_Gbeta/gamma_subunits       1       1
# 8                Ras_activation_uopn_Ca2+_infux_through_NMDA_receptor       1       1
# 9                        Response_to_elevated_platelet_cytosolic_Ca2+       1       1
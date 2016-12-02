univ <- read.csv('caprice/RESULT/11GDS_sam_results.csv', stringsAsFactors = FALSE) %>%
    mutate( gds = paste0('GDS', gds) )

gene2pathwaydf <-
    read.delim("caprice/MAP/UniProt2Reactome.tsv",
               sep='\t', header = TRUE, stringsAsFactors=FALSE) %>%
    mutate( pathway = gsub(' ','_',pathway) )


sampeal <- list() # pea list
peaed_gds <- c()
for( i in seq_along(gdsv) ){ # FIXME gdsv
    message("  ", gdsv[i])
    probev <- univ %>% filter( gds == gdsv[i] ) %>% .$probe
    if (length(probev) < 5){
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
    dplyr::select(pathway, GDS5204, GDS1962) %>% filter( GDS5204 < .001 & GDS1962 < .001 )

#                                                    pathway      GDS5204      GDS1962
# 1  Antigen_processing:_Ubiquitination_&_Proteasome_degradation 1.448784e-06 1.840517e-13
# 2                        Autodegradation_of_Cdh1_by_Cdh1:APC/C 5.051036e-05 5.037400e-04
# 3          Cargo_recognition_for_clathrin-mediated_endocytosis 4.755856e-06 1.287089e-05
# 4         Cdc20:Phospho-APC/C_mediated_degradation_of_Cyclin_A 8.432069e-05 3.130127e-04
# 5                                Clathrin-mediated_endocytosis 1.671509e-08 7.214533e-06
# 6                          COPI-mediated_anterograde_transport 1.986332e-04 7.774106e-06
# 7                              EPHB-mediated_forward_signaling 2.427334e-05 2.417615e-04
# 8        Interactions_of_neurexins_and_neuroligins_at_synapses 8.203021e-04 8.533214e-06
# 9                                               Macroautophagy 2.366191e-04 1.562929e-04
# 10                                       MAPK6/MAPK4_signaling 5.151226e-05 3.865592e-04
# 11                                Mitochondrial_protein_import 1.076748e-05 6.778382e-05
# 12                                    Neutrophil_degranulation 1.614170e-06 1.833281e-07
# 13                                 Orc1_removal_from_chromatin 3.762500e-04 4.132735e-04
# 14             Retrograde_transport_at_the_Trans-Golgi-Network 6.776569e-04 1.966844e-05
# 15                                            Rho_GTPase_cycle 1.010227e-06 2.938500e-05
# 16                             Separation_of_Sister_Chromatids 5.949899e-04 1.801072e-08
# 17                                             UCH_proteinases 8.766066e-05 3.806522e-04

as.data.frame(pea_matrix) %>% rownames_to_column("pathway") %>%
    dplyr::select(pathway, GDS5204, GDS1962) %>% filter( grepl("CaMK|Ca2+",pathway) )
# 1                                          Ca2+_activated_K+_channels 0.49427884 0.0952784214
# 2                                                        Ca2+_pathway 0.01477454 0.0039304532
# 3                            CaMK_IV-mediated_phosphorylation_of_CREB 0.19629016 0.2466968021
# 4               CREB_phosphorylation_through_the_activation_of_CaMKII 0.23303988 0.0009090653
# 5                CREB_phosphorylation_through_the_activation_of_CaMKK 0.07216189 0.1547041681
# 6                                  Elevation_of_cytosolic_Ca2+_levels 0.87529720 0.4138199428
# 7 Inhibition__of_voltage_gated_Ca2+_channels_via_Gbeta/gamma_subunits 0.02940185 0.0051694082
# 8                Ras_activation_uopn_Ca2+_infux_through_NMDA_receptor 0.07947729 0.0003569884
# 9                        Response_to_elevated_platelet_cytosolic_Ca2+ 0.63318160 0.2466968021
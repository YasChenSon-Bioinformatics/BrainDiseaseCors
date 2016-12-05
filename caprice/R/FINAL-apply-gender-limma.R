source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv = c('not-GPL570', 'blacklist'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl, suppress = FALSE)
Ml <- extractMatrixFromEset(ESETl)


topped_nong <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=nrow(Ml[[1]]), p_threshold=.05, method='fdr', design = 'binary')
topped_g    <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene=nrow(Ml[[1]]), p_threshold=.05, method='fdr', design = 'gender')

deg_matrix_g <- build_deg_matrix(topped_g, type = 'F')
deg_matrix_g
# GENDER         GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# (gds ) GDS5204      NA       0       0       0       0       0       0       0       0       0       0
# (gds ) GDS4522       9      NA       0       0       0       0       0       0       0       0       0
# (gds ) GDS4523      30      39      NA       0       0       0       0       0       0       0       0
# (Meta) GDS4358      21       7      19      NA       0       0       0       0       0       0       0
# (    ) GDS4218       7       5       4       0      NA       0       0       0       0       0       0
# (Meta) GDS4136      25       8      29      20       5      NA       0       0       0       0       0
# (    ) GDS4135       1       4       3       2       4       4      NA       0       0       0       0
# (gds ) GDS2821      24       7      32      24       3      30       2      NA       0       0       0
# (    ) GDS2795       8       1       6       1       4       0       5       2      NA       0       0
# (    ) GDS1962      12       3       5       4       2       6       4       4       3      NA       0
# (    ) GDS1917       1       2       1       3       2       5       6       5       5       1      NA

deg_matrix_nong <- build_deg_matrix(topped_nong, type = 't')
deg_matrix_nong
#                GDS5204 GDS4522 GDS4523 GDS4358 GDS4218 GDS4136 GDS4135 GDS2821 GDS2795 GDS1962 GDS1917
# (gds ) GDS5204      NA       0       0       0       0       1       0       0       1       1       0
# (gds ) GDS4522       0      NA      11       4       1       0       1       5       1       3       0
# (gds ) GDS4523       1       5      NA       1       2       0       2       5       2       0       0
# (Meta) GDS4358       0       0       0      NA       0       2       2       2       3       2       0
# (    ) GDS4218       0       5       0       0      NA       0       0       0       1       1       0
# (Meta) GDS4136       5       0       0       0       0      NA       5       1       4       3       3
# (    ) GDS4135       0       0       1       0       0       0      NA       1       5       3       5
# (gds ) GDS2821       4       2       1       0       0       1       0      NA       1       2       1
# (    ) GDS2795       4       0       0       0       0       0       0       1      NA       2       1
# (    ) GDS1962      13       0       2       0       0       3       0       0       1      NA       1
# (    ) GDS1917       0       1       1       0       0       1       0       0       0       0      NA

# Advantages of gender-limma
#
# Advantage 1
i_GDS <- 4
targetGDS <- GDSl[[ i_GDS ]]@header$dataset_id[1]

thisGDScolumns <- GDSl[[ i_GDS ]]@dataTable@columns %>% mutate( sample = as.character(sample) )
after_join_thisGDS <- left_join(getGenderMetaInfo(GDSl[[i_GDS]]@header$dataset_id[1]) %>%
                                    filter( sample %in% thisGDScolumns$sample ),
                                thisGDScolumns, by = 'sample') %>%
    arrange(disease.state) %>% mutate( r = row_number() )

tmp <- topped_g[[i_GDS]]$table %>%  rownames_to_column('probe') %>% 
    #mutate( marker = dz_TRUE * isFemale_TRUE ) %>%
    mutate( marker = abs(isFemale_TRUE) ) %>%
    arrange(desc(marker))

genderProbev <- tmp$probe[1:20]

topped_g[[4]]$table %>% rownames_to_column('probe') %>%
    mutate( genderEffect = abs(isFemale_TRUE) ) %>% arrange(desc(genderEffect)) %>% head(20) %>% .$probe -> genderProbev
    
rownames(topped_nong[[4]]$table)[26:50]

topped_nong[[4]]$table[genderProbev, ]

#Ml[[4]][genderProbev, ] %>%
#, 
gender_effects <-
    Ml[[4]][c('224590_at', '221728_x_at'), ] %>%
    as.data.frame %>% # dplyr cannot handle a class 'matrix'
    rownames_to_column('probe') %>%
    gather(key=smpl, value=eval, -probe) %>% # See: https://blog.rstudio.org/2014/07/22/introducing-tidyr/
    left_join(. , after_join_thisGDS, c('smpl' = 'sample')) %>% # we need disease.state later
    #filter( gender == 'M') %>%
    filter( tissue == 'Frontal cortex' ) %>% mutate( r = row_number() ) %>% mutate( disease = grepl('HIV', disease.state) ) %>%
    group_by(probe,disease,gender) %>% mutate( hypo_Monly = mean(eval) ) %>% 
    group_by(probe,disease) %>% mutate( hypothesis = mean(eval), x = min(r), xend = max(r), hypo_Monly = min(hypo_Monly) )
    

t.test(gender_effects %>% filter( disease & probe == '221728_x_at') %>% .$eval,
       gender_effects %>% filter(!disease & probe == '221728_x_at') %>% .$eval)

png('RESULT/gender-limma-bias.png', width=640, height = 640)
    gender_effects %>% filter(probe == '224590_at') %>%
    ggplot() +
    geom_text(aes(x=r, y=eval, label=gender, color=disease), fontface='bold', size=7) +
    geom_segment(aes(x= x, xend=xend, y=hypothesis,yend=hypothesis, color=disease ), size=1.5) +
    geom_segment(aes(x= x, xend=xend, y=hypo_Monly,yend=hypo_Monly, color=disease), linetype="dashed", size=1.5) +
    #facet_wrap( ~ probe, nrow =5) + # create subplots for each probe
    theme( axis.title.y = element_text(angle=0), text = element_text(face='bold', size=10),
           axis.text.x = element_blank(), axis.text.y = element_text(size=20),
           legend.position = "bottom") +
    labs(title=paste0('Expression Values of Probe "224590_at" in ', targetGDS, ' (HIV study)'),
         subtitle='Horizontal Bars are Hypotheses with Female (solid, p<.01) and without Female (dashed, p=.63)',
         x = 'Sample Index', y = 'Expression\nValues') + scale_y_continuous(breaks=seq(2,3.6,by=.1))
dev.off()

png('RESULT/gender-limma-without-gender-info.png', width=640, height = 640)
gender_effects %>% filter(probe == '224590_at') %>%
    ggplot() +
    geom_point(aes(x=r, y=eval,color=disease), fontface='bold', size=7) +
    geom_segment(aes(x= x, xend=xend, y=hypothesis,yend=hypothesis, color=disease ), size=1.5) +
    theme( axis.title.y = element_text(angle=0), text = element_text(face='bold', size=10),
           axis.text.x = element_blank(), axis.text.y = element_text(size=20),
           legend.position = "bottom") +
    labs(title=paste0('Expression Values of Probe "224590_at" in ', targetGDS, ' (HIV study)'),
         subtitle='Horizontal Bars are Hypotheses (p<.01)',
         x = 'Sample Index', y = 'Expression\nValues') + scale_y_continuous(breaks=seq(2,3.6,by=.1))
dev.off()


topped_nong[[4]]$table['214218_s_at',]
topped_g[[4]]$table['214218_s_at',]


# analysis 1: remove female samples 

all <- nrow(Ml[[1]])
topped_withg   <- applyTtestToGeneExpressionMatrices(Ml,GDSl,nTopGene=all,p_threshold=.05,removeGenderEffect = FALSE)
topped_removeg <- applyTtestToGeneExpressionMatrices(Ml,GDSl,nTopGene=all,p_threshold=.05,removeGenderEffect = TRUE)

gender_results <- 
left_join(
    topped_withg[[4]]$table   %>% rownames_to_column('probe') %>% arrange(P.Value) %>% mutate(rw=row_number()),
    topped_removeg[[4]]$table %>% rownames_to_column('probe') %>% arrange(P.Value) %>% mutate(rr=row_number()),
    by = 'probe'
)

png('caprice/RESULT/GDS4358-1region.png', width=960, height=960)
with(gender_results,
     plot(P.Value.x, P.Value.y, pch=18,cex=.4,
            xlab="P values of 24 samples in GDS4358 (i.e. all samples at Frontal Cortex)"),
            ylab="P values of 23 samples (excluding a female) in GDS4358")
dev.off()

gender_results %>% mutate( px = round(P.Value.x,3), py = round(P.Value.y,3) ) %>% group_by(px,py) %>% summarize( n = n()) %>% ggplot() + geom_point(aes(px,py), alpha=.5)

gender_results %>% dplyr::rename( px = P.Value.x, py = P.Value.y ) %>% filter(py<.05) %>% ggplot() + geom_point(aes(px,py), alpha=.5)

gender_results %>% filter( P.Value.x > .05 & P.Value.y < .05 )
# >   400 probes are newly significant when considering three regions
# > 17900 probes are newly significant when considering one region (Frontal Cortex)

gender_results %>% filter( P.Value.y < .05 ) %>% arrange(adj.P.Val.y)

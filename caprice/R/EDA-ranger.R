library(ranger)

setGlobalConstantList()
loadLibraries()

GDSl <- download_GDSs(skipv=c('not-GPL570', 'no-disease.state', 'no-control'))
ESETl <- convertGDS2ESET(GDSl)
Ml <- extractMatrixFromEset(ESETl)

timev <- c()

i <- 1
vil <- list()
for( i in seq_along(Ml) ){
  univ_unscaled <- t(Ml[[i]])
  univ <- scale(univ_unscaled, center = TRUE, scale = TRUE)  %>% as.data.frame(., stringsAsFactors = FALSE)
  univ$y <- as.factor(GDSl[[i]]@dataTable@columns$disease.state)
  
  start <- proc.time()
  # DO NOT USE FORMULA. protect() overflow error occurred.
  rangered <- ranger(data = univ, num.trees = 10000, mtry = round(sqrt(ncol(univ)),0), min.node.size = 1,
                     scale.permutation.importance = FALSE, importance = "permutation", write.forest = FALSE,
                     save.memory = FALSE, dependent.variable.name = 'y', classification = TRUE, seed=1234)
  end <- proc.time()
  timev <- c(timev,  (end - start)['elapsed'])
  vil[[i]] <- list( vi = rangered$variable.importance, err = rangered$prediction.error, cm = rangered$confusion.matrix)
  message(i)
}

d <- data.frame( probe = names(vil[[2]]$vi), GDS4523 = vil[[2]]$vi, stringsAsFactors = FALSE )
i<-1
for( i in seq_along(vil) ){
  this <- data.frame( probe = names(vil[[i]]$vi), GDS = vil[[i]]$vi, stringsAsFactors = FALSE )
  colnames(this) <- c('probe', GDSl[[i]]@header$dataset_id[1])
  if ( i == 2 )
    next
  d <- left_join(d, this, by='probe')
}
d %>% select(probe, GDS4523, GDS4522) %>% filter( ! is.na(GDS4522) ) %>% filter(GDS4523 > 0 & GDS4522 > 0) %>% mutate(score = GDS4523 * GDS4522) %>% arrange(desc(score)) 

d %>% select(probe, GDS4523) %>% arrange(desc(GDS4523))

function(){
  i <- 2
  GDSl[[2]]@header$dataset_id # GDS4523
  a <- exprs(ESETl[[i]])[ rownames(ESETl[[i]]) %in% c('200011_s_at','203586_s_at', '221482_s_at'), ] %>% as.data.frame(.,stringsAsFactors=FALSE) %>% rownames_to_column('probe') %>% gather(sample,eval, -starts_with('probe')) %>% mutate( sample = factor(sample)) %>% left_join(., GDSl[[i]]@dataTable@columns, by='sample')
  
  t.test( (a %>% filter(disease.state=='control' & probe != '200011_s_at'))$eval,
          (a %>% filter(disease.state!='control' & probe != '200011_s_at' ))$eval )
  
  ggplot(a,aes(x=sample,y=eval,color=disease.state)) +  geom_text(aes(label=sample),size=3,lineheight=.6)+ facet_wrap(  ~  probe, ncol=5)
  
  a %>% group_by(probe, disease.state) %>% summarize( avg = mean(eval) ) %>% 
  
  t.test( )
  
  GDSl[[2]]@header$dataset_id # GDS4523
  GDSl[[2]]@header$sample_count #51
  sss <- GDSl[[2]]@dataTable@table[ GDSl[[2]]@dataTable@table$ID_REF == '200011_s_at', c(-1, -2) ]
  plot( 1:length(sss), y=log2(as.numeric(sss)), col=as.factor(GDSl[[2]]@dataTable@columns$disease.state))
}



function(){
  i <- 1
  candidx <- 1:30
  cands <- names(sort(vil[[i]]$vi, decreasing = TRUE)[candidx])
  tmpdf <- data.frame( probe = cands, imp = sort(vil[[i]]$vi, decreasing = TRUE)[candidx], stringsAsFactors = FALSE )
  material <- exprs(ESETl[[i]])[ rownames(ESETl[[1]]) %in% cands, ] %>% as.data.frame(.,stringsAsFactors=FALSE) %>% rownames_to_column('probe') %>% gather(sample,eval, -starts_with('probe')) %>% mutate( sample = factor(sample)) %>% left_join(., GDSl[[i]]@dataTable@columns, by='sample') %>% mutate( age = gsub('[^0-9]','',age) %>% as.numeric, meta = paste0(age,'',substr(gender,1,1)) ) %>% left_join(., tmpdf, by='probe') %>% mutate( strip = paste0(probe, '  VarImp:', round(imp,9) ) ) 
  ggplot(material,aes(x=sample,y=eval,color=disease.state)) +  geom_text(aes(label=meta),size=3,lineheight=.6) + facet_wrap(  ~  strip, ncol=5) + theme(legend.position="bottom", axis.text.x = element_text(angle=90)) + labs(title = paste0(GDSl[[i]]@header$title,'  -  ', GDSl[[i]]@header$dataset_id[1]), subtitle = strwrap(GDSl[[i]]@header$description[1], width = 140) %>% paste(.,collapse='\n') )
  
  a <- 
  exprs(ESETl[[i]])[ rownames(ESETl[[1]]) %in% c('45288_at','38703_at'), ] %>% as.data.frame(.,stringsAsFactors=FALSE) %>% rownames_to_column('probe') %>% gather(sample,eval, -starts_with('probe')) %>% mutate( sample = factor(sample)) %>% left_join(., GDSl[[i]]@dataTable@columns, by='sample')
  
  t.test( (a %>% filter(disease.state=='control' & probe != '38703_at'))$eval,
          (a %>% filter(disease.state!='control' & probe != '38703_at' ))$eval )
  
  ggplot(a,aes(x=sample,y=eval,color=disease.state)) +  geom_text(aes(label=probe),size=3,lineheight=.6)+ facet_wrap(  ~  probe, ncol=5)
}


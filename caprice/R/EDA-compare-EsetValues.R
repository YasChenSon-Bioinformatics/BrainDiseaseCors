function(){
  expl <- sapply(all_ESET, function(x){ exprs(x) } )
  pdf("/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/EDA/EsetVal1stEd.pdf")
  l <- expl[[1]]
  for( l in expl ) {
    print(hist(l))
  }
  dev.off()
}

function(){
  tmpdf <- sapply(all_GDS, function(x){x@dataTable@table[, 4] %>% as.numeric} )  %>% as.data.frame
  colnames(tmpdf) <- sapply(all_GDS, function(x){x@header$dataset_id[1]} )
  tmpdf <- tmpdf %>% gather(gds,eval,starts_with('GDS'))
  pdf("/Users/PCUser/Dropbox/CU2016/F16CLASSES/TB_Tatonetti/BrainDiseaseCors/caprice/EDA/expValues-logOrNonLog1stEd.pdf")
  for( g in unique(tmpdf$gds) ) {
    p <- ggplot(tmpdf %>% filter(gds == g )) + geom_density(aes(x=eval), fill="red", alpha=.5) + labs(title=g)
    print(p)
  }
  dev.off()
}
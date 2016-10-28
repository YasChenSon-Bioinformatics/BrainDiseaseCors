function(){ # Pre-Analysis: What experimental condition each GDS has?
  # experimental condition vector
  exp_condv <- sapply( all_GDS, function(gds) { colnames(gds@dataTable@columns) }) %>% unlist %>% unique 
  cond_df <- matrix(data = 0, # all-zero amtrix
                    nrow = length(all_GDS),
                    ncol = length(exp_condv),
                    dimnames = list(GDS_on_GPL570v, exp_condv)) %>% data.frame
  for ( i in seq_along(all_GDS) ){
    exp_cond_of_this_dataset <- colnames(all_GDS[[i]]@dataTable@columns)
    cond_df[i, ] <- as.numeric(exp_condv %in% exp_cond_of_this_dataset)
  }
  cdf <- cond_df %>% add_rownames(var='GDS') %>% gather(key = meta, value = bit, -starts_with('GDS'))
  experimental_bitmap <-
    left_join(cond_df %>% add_rownames(var='GDS') %>% gather(key = meta, value = bit, -matches('GDS')),
              cdf %>% group_by(meta) %>% summarize( nentry = sum(bit) ), by="meta") %>%
    mutate(meta=paste(nentry,meta)) %>% ggplot(aes(x=gsub('GDS','',GDS), y=reorder(meta,nentry))) +
    geom_tile(aes(fill=as.factor(bit))) + scale_fill_grey()
  print(experimental_bitmap)
}

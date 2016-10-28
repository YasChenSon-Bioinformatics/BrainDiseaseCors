common_probe_ids = rownames(all_M[[1]])
for (m in all_M) {
   common_probe_ids <- intersect(common_probe_ids, rownames(m))
   message(length(common_probe_ids))
}
affyIDs <- lapply( all_GDS, function(x) { x@dataTable@table$ID_REF} ) %>% unlist %>% unique
clip <- pipe("pbcopy", "w")
write.table(paste(affyIDs, collapse='\n'), file=clip)
close(clip)
# NB: 54677 Common Columns when applied as.data.frame,
#   which contains 'sample' and 'description' as last 2 columns
# When applied to matrix, 54675 probe ids are common.
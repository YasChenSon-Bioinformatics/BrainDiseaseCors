source('/home/PCUser/BrainDiseaseCors/caprice/R/analyze-GPL570.R')

setGlobalConstantList(rootDir = '/home/PCUser/BrainDiseaseCors')
loadLibraries()
GDSl <- donwload_GDSs()

# Note: some rows like 1320_at or 1405_i_at in GDS4522 are <NA>.
# I asked The GEO team about this issue. And <NA> is for unreliable expression values.
GDS_on_GPL570v <- sapply( all_GDS, function(gds) {gds@header$dataset_id[1]} )

ESETl <- convertGDS2ESET(GDSl)


sapply(all_M, dim)
sapply(all_M, function(x) {sum(is.na(x))} )  # MAYBE-LATER What's these nas?
sapply(all_GDS, function(x) {message(paste0(unique(x@dataTable@columns$disease.state),collapse="\t"))} )  # MAYBE-LATER What's these nas?

Ml <- extractMatrixFromEset(ESETl)
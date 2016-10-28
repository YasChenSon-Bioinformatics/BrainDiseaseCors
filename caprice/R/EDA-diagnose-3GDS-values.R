# Limma provides lots of diagnostics functions.

# First we need to understand various Biobase/limma data structures we can apply diagnostics to:
#
#  'ExpressionSet' is in the Biobase package. Its purpose is the same as 'EList' below.
#  see: help('ExpressionSet-class')
#
#  'EList'  : a list of normalized expression values (E-values).
#             one-channel microarrays.
#             'EListRaw' has expression values on the raw scale from read.maimages().
#             'EList' has expression values on the log scale,
#             usually after background correction and normalization by normalizeBetweenArrays() or by voom().
#             For further info: help('EList-class').
#
#   FYI
# 
#  'RGList' : list of red/green channel foreground/background intensities. help('RGList-class')
#  'MAList' : list of M-values (log2 expression ratios) and A-values (average log2 expression values). help('MAList-class')
#

# "Plots for individual arrays include the foreground–background plots mentioned above
#  (plotFG), image plots that can reveal inconsistencies across the array surface
#  (imageplot) and mean-difference plots that show intensity-dependent trends in the 
# log-ratios of two-colour arrays (plotMA, Figure Figure3B).3B). The plotMA function 
# can show similar plots for single channel data."
# (from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/)

# Tool 1: plotMA(). Useful for what? I'll check later
# MA plot is just a kind of scatterplot. But both x-axis and y-axis need to be computed from the data.
#
# MA plot is first proposed as a method to investigate two-color microarray data,
# but it can be actually used for single-channel data too
#
# M >  1 means 'expression level increased' (logFC >  1).
# M == 0 means 'expression level is not changed' (logFC == 0).
# M < -1 means 'expression level decreased' (logFC < -1).
#
plotMA( ESETl[[1]], array=1, xlab="Average log-expression (A-values)", ylab="Expression log-ratio (this sample vs others, M-values)")
plotMA( ESETl[[1]], array=2 )
plotMA( ESETl[[2]] ) # MAYBE-LATER: why linearly correlated?

# Tool 2: plotDensities(). Useful for checking some samples have irregular value ranges or not.
# x-axis: exp level
# y-axis: relative frequency
# color: sample
plotDensities( ESETl[[1]], legend = c("topleft", "topright")[1] )

# What does plotDensities() do? Reproduction by R default functions is like
hist( exprs(ESETl[[1]])[, 'GSM439800'], breaks =100, prob=TRUE, border="white" ) # draw plot area
for( l in list( first = c('GSM439800', 'red'), second  = c('GSM439790', 'blue')) ){
  lines( density( exprs(ESETl[[1]])[ , l[1] ] ), col=l[2])
  # hist( exprs(ESETl[[1]])[, 'GSM439800'], breaks =100, prob=TRUE ) # exp. vals for GSM439800
}

# Other experiments

plot(exprs(ESETl[[1]])[,1], exprs(ESETl[[1]])[,2], pch=18)
# Summary: gene expression comparison between 2 samples.
# x-axis: exp levels of sample 1
# y-axis: exp levels of sample 2
# Interpretation: mostly linearly correlated (since they measure the same genes).
#                 But if sample 1 is healthy while 2 is diseased and if points digress from linear relationship, 
#                 the gene might be related to the disease.

boxplot( exprs(ESETl[[1]]), pch=18)
boxplot( exprs(ESETl[[2]]), pch=18) # strange.
# Summary: gene expression range comparison for all samples.
# x-axis: sample name
# y-axis: exp level
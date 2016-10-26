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
options(mfrow=c(2,2))
plotMA( ESETl[[1]], array=1, xlab="Average log-expression (A-values)", ylab="Expression log-ratio (this sample vs others, M-values)")
plotMA( ESETl[[1]], array=2 )
plotMA( ESETl[[2]] ) # why linearly correlated?

# Tool 2: plotDensities(). Useful for checking some samples have irregular value ranges or not.
plotDensities( ESETl[[1]], legend = c("topleft", "topright")[1] )

# What does plotDensities() do? Reproduction by R default functions is like
hist( exprs(ESETl[[1]])[, 'GSM439800'], breaks =100, prob=TRUE, border="white" ) # draw plot area
for( l in list( first = c('GSM439800', 'red'), second  = c('GSM439790', 'blue')) ){
  lines( density( exprs(ESETl[[1]])[ , l[1] ] ), col=l[2])
  # hist( exprs(ESETl[[1]])[, 'GSM439800'], breaks =100, prob=TRUE ) # exp. vals for GSM439800
}


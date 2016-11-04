# setup

source('BrainDiseaseCors/caprice/R/analyze-GPL570.R')
setGlobalConstantList(rootDir = '/Users/admin/Dropbox/Columbia2016/Bioinformatic/Projects/Project_SourceTree')
loadLibraries()

# search Differentially Expressed Genes for each datasets 

 GDSl <- download_GDSs(GDSnumberv = c(4523, 4522, 4358, 4218, 4136, 2821, 1917), skipv = c('not-GPL570'), suppress = FALSE)
ESETl <- convertGDS2ESET(GDSl)
   Ml <- extractMatrixFromEset(ESETl)
topped <- applyTtestToGeneExpressionMatrices(Ml, GDSl, nTopGene = 400, method = 'fdr') # CAUTION: currently our design matrices are just 'disease.state'. no age or gender
DEG_matrix <- sapply(topped, function(x) rownames(x$table))
length(       as.vector(DEG_matrix) ) # 2800
length(unique(as.vector(DEG_matrix))) # 2720

commonDEGv <- names(table(DEG_matrix)[table(DEG_matrix) > 1])
GDSnamev   <- sapply(GDSl, function(gds){gds@header$dataset_id[1]})

i <- 2
disease_signature_matrix <- matrix(NA, nrow = length(commonDEGv), ncol = length(ESETl), dimnames =list(commonDEGv, GDSnamev) )
rownames(disease_signature_matrix)
for ( i in seq_along(ESETl) ){
  this <- exprs(ESETl[[i]])[commonDEGv, ]
  control_samplev <- grepl('control',GDSl[[i]]@dataTable@columns$disease.state)
  disease_signature_matrix[ , i ] <- rowMeans(this[, control_samplev]) - mean(this[, ! control_samplev], na.rm=TRUE)
  # NB: GDS4522, GDS4218, GDS2821 is not on GPL570 so na.rm=FALSE is necessary
}
disease_signature_matrix <- as.matrix(na.omit(as.data.frame(disease_signature_matrix))) # remove rows with na 

boxplot(disease_signature_matrix) # FIXME: seems that we need to normalize value ranges

###########################################################
# Compute correlation between diseases expression values  #
###########################################################

###########################################################
# Sd normal  #
###########################################################
column_sd <- apply(disease_signature_matrix, 2, sd)
column_mean <- apply(disease_signature_matrix, 2, mean)
normalized_disease_signature_matrix_new <- (disease_signature_matrix - column_mean)/column_sd


###########################################################
# Min-Max normal  #
###########################################################
standardize <- function(z) {
  rowmax <- apply(z, 2, max) #max
  rowmin <- apply(z, 2, min)  # min absolute deviation
  denomenator <- rowmax - rowmin
  nomenator <- sweep(z, 2, rowmin,"-")  #subtracting min expression
  result <- sweep(nomenator, 2, denomenator, "/")  # dividing by median absolute deviation
  return(result)
}

normalized_matrix <- standardize(disease_signature_matrix)


if(require('corrplot') == FALSE) install.packages("corrplot"); library(corrplot)

#cv = cor(normalized, use="complete.obs") # Cor vector
cv = cor(normalized_matrix, use="complete.obs")
namev = c("GDS4523 - \nSchizophrenia",
          "GDS4522 - \nSchizophrenia",
          "GDS4358 - \nHIV",
          "GDS4218 - \nMultiple \nSclerosis",
          "GDS4136 - \nAlzheimers",
          "GDS2821 - \nParkinsons",
          "GDS1917 - \nSchizophrenia")
colnames(cv) <- rownames(cv) <- namev

#attach(cv)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
#corrplot.mixed(cv, t1.pos="r", t1.col="blue", c1.srt=60, c1.pos="r", cl.align.text="r", mar=c(1,1,1,1), height=1600, width=1600)
# perhaps you'd like to do this?
corrplot.mixed(cv, tl.pos=c("d", "lt", "n")[1], tl.col="blue", cl.align.text="r", mar=c(1,1,1,1))

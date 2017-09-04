# Brain Disease Correlation Analysis Project

+ Fall 2016 Bioinformatics class project
+ [Final Report](https://github.com/YasChenSon-Bioinformatics/BrainDiseaseCors/blob/master/writeups/Final_Writeup.pdf)

# RESOURCES

+ [2-page proposal](https://docs.google.com/document/d/1WH9bjXNLgi4JiFfaLSqGhYR2SLK-xyDZ1bOGP8bEDcI/edit)
+ [Microarray data analysis introduction (PH525x)](http://genomicsclass.github.io/book/)


# Facts worthy of attention
+ **Some incorrect expression values in GSE are transformed into null in GDS**. Although GDS (Datasets) uses GSE (Series) as data sources, some expression values are different from GSE. In other words GDS is *"curated"* by The GEO team. So GDS is much helpful for determining which rows to be used.
+ **Some GDS needs do.log2=TRUE, and others need do.log2=FALSE.** NaN in your expression values might be due to multiple log2 transformations.
+ **getGEO() might return partial .tar.gz files without any warnings.** Using locally-cached files is strongly recommended.
+ **Some probe IDs of GPL570 and GPL96 are common.** Comparison might be possible.
+ **There seems to be no 1-to-1 maps between Affymetrix Probe IDs and various types of gene IDs**. I tried Entrez Gene ID, Ensembl ID, and Gene symbols. All failed. The better way might be to use Affymetrix Probe IDs as is.

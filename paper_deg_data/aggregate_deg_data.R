"
Module for confirming our results from the SAM experiments with the data from the papers
We have only organized paper reusults from four of the datasets so far
"

require(gdata)
setwd('~/Desktop/Columbia/Bioinformatics/project/paper_deg_data/')

datasets_num = list(1917, 1962, 2795, 2821, 4135, 4136, 4218, 4358, 4522, 4523, 5204)
datasets_w_xls = list(2821, 4136, 4218, 4358, 4522, 4523) # include 2431?

colnames_2821 = c("Gene Symbol", "Gene Name", "FragmentID(ChipID)", "Probe Set ID", "t-test p-Value sn", "Fold change sn", "% Present Cases sn", "% Present controls sn", "t-test p-Value putamen", "Fold change putamen", "% Present Cases putamen", "% Present controls putamen", "t-test p-Value caudate", "Fold change caudate", "% Present Cases caudate", "% Present controls caudate")
colnames_4136 = c("Probe set", "Symbol", "Description", "r overall MMSE", "p overall MMSE", "r overall NFT", "p overall NFT", "r incip. MMSE", "p incip. MMSE", "r incip NFT", "p incip NFT")
colnames_4218 = c("Probe Set ID", "Representative Public ID", "Gene title", "Control 1", "Control 2", "CP1", "CP2", "CAP1", "CAP2", "AP")
colnames_4231 = c("Probe Set ID up", "Unigene up", "Gene Title up", "Gene Symbol up", "FC up", "pv up", "empty", "Probe Set ID down", "Unigene down", "Gene Title down", "Gene Symbol down", "FC down", "pv down")
colnames_4358_first_5_sheets = c("row#names", "GeneName", "GeneIndex", "foldChange", "AntiLogFC", "Norm FC", "Pvalue", "testStat", "AdjPvalue", "Signif#p", "Gene Symbol", "Gene Title", "Pathway", "go biological process term", "go molecular function term", "go cellular component term")  
colnames_4358_last_4_sheets = c("Affymetrix #", "Gene Symbol", "Gene Title", "Affymetrix #", "Fold Change", "P value", "BH adjusted P", "GeneIndex", "AntiLog FC", "Norm FC", "testStat", "Signif#p", "Pathway", "go biological process term", "go molecular function term", "go cellular component term")
colnames_4523_first_sheet = c("Affy ID", "HUGO sumbol", "EntrezGene ID", "Gene Description", "HBB Model Fold Change", "HBB Model P Value", "HBB QC Covariate PV", "HBB Model Age PV", "HBB Model Disease PV", "HBB Model Gender PV", "HBB Control Median Intensity", "HBB Schizo Median Intensity")
colnames_4523_second_sheet = c("Affy ID", "Hugo Symbol", "EntrezGene ID", "Gene Description", "CCHPC Model Fold Change", "CCHPC Model P Value", "Control Median Intensity (Raw)", "Schizo Median Intensity (Raw)", "CCHPC Model Gender PV", "CCHPC Model Disease*Gender PV", "CCHPC Model Age PV", "CCHPC Model Age*Gender PV", "CCHPC QC Covariate PV")
colnames_4523_third_sheet = c("Affy ID", "Gene Name", "HUGO Symbol", "Entrez Gene ID", "Gene Description", "CCHPC Fold Change", "CCHPC P value", "HBB Fold Change", "HBB P Value")

# Todo - figure out how to access the sheets from the xls files
papers_to_confirm = list(2821, 4136, 4218)
colnames_list = list(colnames_2821, colnames_4136, colnames_4218)
gene_name_columns = c('') # ... todo ...

# Need:
# 1917 - paper not found
# 1962 - no supplemental data (SCF is a DEG though)
# 2795 - supp not included - email authors?
# 4135 - .doc files --> how to export the tables to xls?
# 5204 - no supp data - email authors?

i = 1
papers_deg_data = list()
for (num in papers_to_confirm) {
  message(i, num, ' --------------------------')
  papers_deg_data[[i]] = read.xls(paste(num, '_supp.xls', sep=''), header=TRUE, blank.lines.skip=TRUE)
  colnames(papers_deg_data[[i]]) = colnames_list[[i]]
  i = i + 1
}

# Find intersections (todo - make indices match):
common_degs = list()
# Add "_at" to our deg's so we can compare to paper degs
# "35179_at" "34221_at" "52731_at" "51192_at" "35626_at"
i = 1
for (df in papers_deg_data) {
  common_degs[[i]] = intersect(lapply(deg_up[[i]], paste0, "_at"), papers_deg_data[[i]]$X)
  i = i + 1
}


intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)
intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)
intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)

# Checking that all are degs in paper deg data 2821
tmp = papers_deg_data[[1]]
tmp[which(tmp$X.4 < .05 | tmp$X.8 < .05 | tmp$X.10 < .05),]
tmp$X.4 = as.numeric(tmp$X.4)
tmp$X.8 = as.numeric(tmp$X.8)
tmp$X.10 = as.numeric(tmp$X.10)

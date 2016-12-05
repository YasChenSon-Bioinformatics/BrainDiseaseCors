"
Module for confirming our results from the SAM experiments with the data from the papers.
Some paper results are still missing or unused

# Still Need:
# 1917 - paper not found
# 1962 - no supplemental data (SCF is a DEG though)
# 2795 - supp not included - email authors?
# 4135 - .doc files --> how to export the tables to xls?
# 5204 - no supp data - email authors?
"

require(gdata)
setwd('~/Desktop/Columbia/Bioinformatics/project/paper_deg_data/')

datasets_num = list(1917, 1962, 2795, 2821, 4135, 4136, 4218, 4358, 4522, 4523, 5204)
datasets_w_xls = list(2821, 4136, 4218, 4358, 4523) # include 2431? 4522?

colnames_2821 = c("Gene Symbol", "Gene Name", "FragmentID(ChipID)", "Probe id", "t-test p-Value sn", "Fold change sn", "% Present Cases sn", "% Present controls sn", "t-test p-Value putamen", "Fold change putamen", "% Present Cases putamen", "% Present controls putamen", "t-test p-Value caudate", "Fold change caudate", "% Present Cases caudate", "% Present controls caudate")
colnames_4136 = c("Probe id", "Symbol", "Description", "r overall MMSE", "p overall MMSE", "r overall NFT", "p overall NFT", "r incip. MMSE", "p incip. MMSE", "r incip NFT", "p incip NFT")
colnames_4218 = c("Probe id", "Representative Public ID", "Gene title", "Control 1", "Control 2", "CP1", "CP2", "CAP1", "CAP2", "AP")
colnames_4231 = c("Probe Set ID up", "Unigene up", "Gene Title up", "Gene Symbol up", "FC up", "pv up", "empty", "Probe Set ID down", "Unigene down", "Gene Title down", "Gene Symbol down", "FC down", "pv down")

colnames_4358_first_5_sheets = c("row#names", "GeneName", "GeneIndex", "foldChange", "AntiLogFC", "Norm FC", "Pvalue", "testStat", "AdjPvalue", "Signif#p", "Gene Symbol", "Gene Title", "Pathway", "go biological process term", "go molecular function term", "go cellular component term")  
colnames_4358_last_4_sheets = c("Probe id", "Gene Symbol", "Gene Title", "Affymetrix #", "Fold Change", "P value", "BH adjusted P", "GeneIndex", "AntiLog FC", "Norm FC", "testStat", "Signif#p", "Pathway", "go biological process term", "go molecular function term", "go cellular component term")

colnames_4523_first_sheet = c("Probe id", "HUGO sumbol", "EntrezGene ID", "Gene Description", "HBB Model Fold Change", "HBB Model P Value", "HBB QC Covariate PV", "HBB Model Age PV", "HBB Model Disease PV", "HBB Model Gender PV", "HBB Control Median Intensity", "HBB Schizo Median Intensity")
colnames_4523_second_sheet = c("Probe id", "Hugo Symbol", "EntrezGene ID", "Gene Description", "CCHPC Model Fold Change", "CCHPC Model P Value", "Control Median Intensity (Raw)", "Schizo Median Intensity (Raw)", "CCHPC Model Gender PV", "CCHPC Model Disease*Gender PV", "CCHPC Model Age PV", "CCHPC Model Age*Gender PV", "CCHPC QC Covariate PV")
colnames_4523_third_sheet = c("Probe id", "Gene Name", "HUGO Symbol", "Entrez Gene ID", "Gene Description", "CCHPC Fold Change", "CCHPC P value", "HBB Fold Change", "HBB P Value")

# 6, 5, 4
papers_to_confirm = list(2821, 4136, 4218) # 4522 too large for now
colnames_list = list(colnames_2821, colnames_4136, colnames_4218)
rownames(deg_matrix) = datasets_w_xls

# Our matrix of papers deg counts, sam results deg counts, and overlap count
deg_matrix <- matrix( NA, ncol = 3, nrow = length(datasets_w_xls) )
colnames(deg_matrix) = c('Paper DEGs count', 'SAM results DEG count', 'Overlap')


### Do first four datasets excel files
i = 1
papers_degs = list()
for (num in papers_to_confirm) {
  message(i, ', ', num, ' --------------------------')
  papers_deg_data = read.xls(paste(num, '_supp.xls', sep=''), header=TRUE, blank.lines.skip=TRUE)
  colnames(papers_deg_data) = colnames_list[[i]]
  papers_degs[[i]] = papers_deg_data$`Probe id`
  i = i + 1
}

# Find intersections:
common_degs = list()
# of 2821
common_degs[[1]] = intersect(union(deg_up[[6]], deg_lo[[6]]), papers_degs[[1]])
deg_matrix[1,1] <- length(papers_degs[[1]])
deg_matrix[1,2] <- length(deg_up[[6]]) + length(deg_lo[[6]])
deg_matrix[1,3] <- length(common_degs[[1]])
# of 4136
common_degs[[2]] = intersect(union(deg_up[[5]], deg_lo[[5]]), papers_degs[[2]])
deg_matrix[2,1] <- length(papers_degs[[2]])
deg_matrix[2,2] <- length(deg_up[[5]]) + length(deg_lo[[5]])
deg_matrix[2,3] <- length(common_degs[[2]])
# of 4218
common_degs[[3]] = intersect(union(deg_up[[4]], deg_lo[[4]]), papers_degs[[3]])
deg_matrix[3,1] <- length(papers_degs[[3]])
deg_matrix[3,2] <- length(deg_up[[4]]) + length(deg_lo[[4]])
deg_matrix[3,3] <- length(common_degs[[3]])

### Handle gds 4358 with 9 sheets
paper_degs_4358 = list()
for (i in 1:5) {
  message(i, num, ' --------------------------')
  xls_data = read.xls('4358_supp.xls', header=TRUE, blank.lines.skip=TRUE, sheet=i)
  colnames(xls_data) = colnames_4358_first_5_sheets
  paper_degs_4358[[i]] = xls_data$`Probe id`
}
for (i in 6:9) {
  message(i, num, ' --------------------------')
  xls_data = read.xls('4358_supp.xls', header=TRUE, blank.lines.skip=TRUE, sheet=i)
  colnames(xls_data) = colnames_4358_last_4_sheets
  paper_degs_4358[[i]] = xls_data$`Probe id`
}
union_deg_4358 = list()
for (i in 1:9) {
  union_deg_4358 = union(union_deg_4358, paper_degs_4358[[i]])
}

# Find intersections:
deg_matrix[4,1] <- length(union_deg_4358)
deg_matrix[4,2] <- length(deg_up[[3]]) + length(deg_lo[[3]]) # 4358 is stored in deg_up/lo [[3]]
deg_matrix[4,3] <- length(intersect(union_deg_4358, union(deg_up[[3]], deg_lo[[3]])))


### Handle gds 4523 with 3 sheets
paper_degs_4523 = list()
# First sheet
xls_data = read.xls('4358_supp.xls', header=TRUE, blank.lines.skip=TRUE, sheet=1)
colnames(xls_data) = colnames_4523_first_sheet
paper_degs_4523[[1]] = xls_data$`Probe id`
# Second sheet
xls_data = read.xls('4358_supp.xls', header=TRUE, blank.lines.skip=TRUE, sheet=2)
colnames(xls_data) = colnames_4523_second_sheet
paper_degs_4523[[2]] = xls_data$`Probe id`
# Third sheet
xls_data = read.xls('4358_supp.xls', header=TRUE, blank.lines.skip=TRUE, sheet=3)
colnames(xls_data) = colnames_4523_third_sheet
paper_degs_4523[[3]] = xls_data$`Probe id`
paper_degs_4523_all = union(paper_degs_4523[[1]], list(paper_degs_4523[[2]], paper_degs_4523[[3]]))

# Find intersections:
deg_matrix[5,1] <- length(paper_degs_4523_all)
deg_matrix[5,2] <- length(deg_up[[1]]) + length(deg_lo[[1]]) # 4523 is stored in deg_up/lo [[1]]
deg_matrix[5,3] <- length(intersect(paper_degs_4523_all, union(deg_up[[1]], deg_lo[[1]])))


# Results:
# Paper           DEGs count SAM results DEG count Overlap
# 2821              151                    53       0
# 4136             3466                    33      22
# 4218             6604                     0       0
# 4358             1098                   850     104
# 4523               86                     2       0


write.table(deg_matrix)


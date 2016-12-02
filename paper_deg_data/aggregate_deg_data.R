require(gdata)
setwd('~/Desktop/Columbia/Bioinformatics/project/paper_deg_data/')

datasets_num = list(1917, 1962, 2795, 2821, 4135, 4136, 4218, 4358, 4522, 4523, 5204)

datasets_w_xls = list(2821, 4136, 4218, 4358, 4522, 4523) # include 2431?

# Need:
# 1917 - paper not found
# 1962 - no supplemental data (SCF is a DEG though)
# 2795 - supp not included - email authors?
# 4135 - .doc files --> how to export the tables to xls?
# 5204 - no supp data - email authors?

i = 1
papers_deg_data = list()
for (num in datasets_w_xls) {
  message(i, num, ' --------------------------')
  if (i == 5) { i = i + 1; next } # 5th dataset takes too long...
  papers_deg_data[[i]] = read.xls(paste(num, '_supp.xls', sep=''), header=TRUE)
  i = i + 1
}

# Find intersections (todo - make indices match):
common_degs = list()
# Add "_at" to our deg's so we can compare to paper degs
intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)
# "35179_at" "34221_at" "52731_at" "51192_at" "35626_at"
i = 1
for (df in papers_deg_data) {
  common_degs[[i]] = 
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

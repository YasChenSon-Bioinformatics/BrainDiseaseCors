require(gdata)
setwd('~/Desktop/Columbia/Bioinformatics/project/paper_deg_data/')

datasets_num = list(1917, 1962, 2795, 2821, 4135, 4136, 4218, 4358, 4522, 4523, 5204)

datasets_w_xls = list(2821, 4136, 4218, 4358, 4522, 4523)

# Need:
# 1917 - paper not found
# 1962 - no supplemental data (SCF is a DEG though)
# 2795 - supp not included - email authors?
# 4135 - .doc files --> how to export the tables to xls?
# 5204 - no supp data - email authors?

i = 1
papers_deg_data = list()
for (num in datasets_num) {
  message(i, num, ' --------------------------')
  papers_deg_data[[i]] = read.xls(paste(datasets_w_xls[[i]], '_supp.xls', sep=''))
  i = i + 1
}

# Find intersections (todo - make indices match):
common_degs = list()
# Add "_at" to our deg's so we can compare to paper degs
intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)
# "35179_at" "34221_at" "52731_at" "51192_at" "35626_at"


intersect(lapply(deg_up[[5]], paste0, "_at"), papers_deg_data[[2]]$X)


supp_4522 = read.csv('4522_supp')


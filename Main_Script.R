all_GDSv <- c(5204,4879,4859,4838,4758,4532,4522,4523,4477,4414,4358,4231,4218,4154,4136,4135,3834,3502,3459,3345,3129,3128,3113,3110,3069,2978,2941,2821,2795,2613,2191,2190,2154,1962,1917,1912,1835,1816,1815,1813,1726,1253,1096,1085,910,909,833,810,707,596,564,426,232,181)

pfc = c(2190, 3502, 4414, 4523) #, 4532) - not in same region



# Download all the datasets. Store only datasets on GPL570
all_GDS = list()
j = 1
for (i in all_GDSv){
<<<<<<< HEAD
  all_GDS[[j]] = assign(paste('GDS', all_GDSv[i], sep='_'), getGEO(GEO=paste("GDS", i, sep=''), destdir='/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project/data'))
=======
>>>>>>> ae9e79fcd6b937d94da196434cc28ab466d0baa5
  # assign('x', 1) stores an integer value 1 to a variable called 'x'.
  # y <- assign('x', 1) stores an integer value 1 to both a variable called 'x' and 'y'.
  # Since you used only all_GDS in this R file, no need to use assign().
  if (Sys.info()['user'] == 'PCUser') {
      tmp = getGEO(GEO=paste("GDS", i, sep=''), destdir='/Users/PCUser/Downloads/Rtmp')
  } else {
      tmp = getGEO(GEO=paste("GDS", i, sep=''), destdir='/Users/ianjohnson/Desktop/Columbia/Bioinformatics/project')
  }
  if ( tmp@header$platform != 'GPL570' ){
      message('-----GDS ', i, ' is skipped because the platform is not GPL570')
      next # go to next for-loop
  }
  all_GDS[[j]] = tmp
  j = j + 1
}

# Extract all to ESETs
all_ESET = list()
j = 1
for (gds in all_GDS){
  all_ESET[[j]] <- GDS2eSet(GDS=gds, do.log2 = TRUE)

  if( any(is.nan(exprs(all_ESET[[j]]))) ){
    # Note: log2(x) is undefined and return NaN if x is negative.
    # What's worse, some datasets looks like already applied log2 transformaiton.
    message("    ", gds@header$dataset_id[1], " seems already applied log2(). SKIP")
    all_ESET[[j]] <- GDS2eSet(GDS=gds, do.log2 = FALSE)
  }
  j = j + 1
}
# check the expression value range by this oneliner.
# sapply(all_ESET, function(x) { summary(as.vector(exprs(x))) })

# To df's
j = 1;
for (eset in all_ESET) {
  gdsNo <- gsub('GDS','', eset@experimentData@other$dataset_id[1])
  assign(paste('df', gdsNo, sep='_'), as.data.frame(eset), envir = .GlobalEnv)
  j = j + 1
}


######################################################################
# Notes:

# Certain GDS's have        Funny column names
# 1085, 1813, 3834      "x1", "x2", etc --> what probes are these?
# 1835                  x882, x880
# 3113                  x156427, x139282
# 4758, 4879            X7896736
# 833                   175019765, 175019766
# 909, 910              x559, x558 




# Limit our analyses to GPL570 datasets (most common, and most comprehensive array)
gpl570 = list(df_5204,df_4838,df_4532,df_4523,df_4522,df_4477,df_4358,df_4231,df_4218,df_4136,df_4135,df_3502,df_2821,df_2795,df_2154,df_1962,df_1917)
<<<<<<< HEAD

gpl570_disease = list(df_5204,df_4838,df_4532,df_4523,df_4522,df_4477,df_4358,df_4231,df_4218,df_4136,df_4135,df_3502,df_2821,df_2795,df_2154,df_1962,df_1917)

=======
gpl570_disease = list(df_5204,df_4838,df_4532,df_4523,df_4522,df_4477,df_4358,df_4231,df_4218,df_4136,df_4135,df_3502,df_2821,df_2795,df_2154,df_1962,df_1917)
>>>>>>> ae9e79fcd6b937d94da196434cc28ab466d0baa5
gpl570_with_diesase_state_column = list(df_4838,df_4523,df_4522,df_4358,df_4231,df_4218,df_4136,df_4135,df_3502,df_2821,df_2795,df_2154,df_1962,df_1917)
# Dropped datasets with no "control" in disease.state column
gpl570_controlled_subset = list(df_4523,df_4522,df_4358,df_4231,df_4218,df_4136,df_3502,df_2821,df_1917)
# Dropped two data sets whose control samples did not show sufficient correlation
gpl570_controlled_correlated_subset = list(df_4523,df_4522,df_4358,df_4218,df_4136,df_2821,df_1917)

# Find NAs (NANs are a todo)
nan_count <-function (x) sapply(x, function(y) sum(is.nan(y)))
for (df in gpl570_disease) {
  print(any(is.na(df)))
  # print(any(is.nan(df)))
}


#####################################################################
# Notes:
# df 1, 3, 6 have no disease.state column!
# df_4838 (neuroectodermal tumors) has no control samples

# gpl570 dataset        variable count        NAs?        NANs?
# df_5204                   54679               n           
# df_4838                   54678               n
# df_4532                   54680               n
# df_4523                   45680               n
# df_4522                   54680               y
# df_4477                   54678               y
# df_4358                   54679               n
# df_4231                   54679               y
# df_4218                   54680               y
# df_4136                   54679               n
# df_4135                   54679               y
# df_3502                   54681               y
# df_2821                   54679               y
# df_2795                   54679               n
# df_2154                   54678               n
# df_1962                   54679               y
# df_1917                   54678               n

# Find the probes common to all df's
common_cols = colnames(df_1917)
for (df in gpl570_disease) {
  common_cols <- intersect(common_cols, colnames(df))
}

# --> 54677 Common Columns!



library(stringr)
library(dplyr)

# Make cleaned subsets for control and disease states for each dataset
control_dfs = list()
disease_dfs = list()
i = 1
for (df in gpl570_controlled_correlated_subset) {
  message(i, "-----------------------------------------------")
  str(df)
  # Subset for control instances
  df_disease = df[!grepl("control", df$disease.state),]
  # Subset for DZ instances
  df_control = df[grepl("control", df$disease.state),]
  
  # todo - are these next two sections working?
  # Remove uncommon columns
  df_control = df_control[, (names(df_control) %in% common_cols)]
  df_disease = df_disease[, (names(df_disease) %in% common_cols)]
  # Remove metadata
  drops <- c("sample", "age", "gender", "tissue", "genotype.variation", "development.stage", "agent", "other", "cell.type", "disease.state", "description", "individual")
  df_control = df_control[, !(names(df_control) %in% drops)]
  df_disease = df_disease[, !(names(df_disease) %in% drops)]
  # Save cleaned control and disease data frames in the list
  
  control_dfs[[i]] = df_control
  disease_dfs[[i]] = df_disease
  i = i + 1
}


#######################################################################
# Note the dimensions of the data frames (separated by ctrl/dz)
##################### Controls
# [1]    23 54680
# [1]    19 54680
# [1]    18 54679
# [1]     9 54679
# [1]     2 54680
# [1]     8 54679
# [1]     0 54679   no control
# [1]    11 54681
# [1]     9 54679
# [1]     0 54679   no control
# [1]     0 54678   no control
# [1]     0 54679   no control
# [1]    14 54678
###################### DZ
# [1]    28 54680
# [1]    23 54680
# [1]    54 54679
# [1]    26 54679
# [1]     5 54680
# [1]    22 54679
# [1]    18 54679   no control
# [1]    21 54681
# [1]    16 54679
# [1]    20 54679   no control
# [1]    12 54678   no control
# [1]   180 54679   no control
# [1]    14 54678
###################### Total
# [1]    51 54680
# [1]    42 54680
# [1]    72 54679
# [1]    35 54679
# [1]     7 54680
# [1]    30 54679
# [1]    18 54679   no control
# [1]    32 54681
# [1]    25 54679
# [1]    20 54679   no control
# [1]    12 54678   no control
# [1]   180 54679   no control
# [1]    28 54678
##################



# Compute correlations within control samples from each dataset
# todo - what happened to control_dfs [0]? --> no control samples!
cors = list()
i = 1
for (df in control_dfs) {
  print(i)
  print(dim(df))
  print("-----------------------------------------------")
  str(df)
  # Transpose to use cor()
  df_t = t(df)
  # Compute corr
  corre = cor(df_t, use="complete.obs")
  cors[[i]] = corre
  i = i + 1
}
names(cors) = c('df_4523','df_4522','df_4358','df_4231','df_4218','df_4136','df_3502','df_2821','df_1917')
attach(cors)
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
plot(cors[[1]], main="df_4523")
plot(cors[[2]], main="df_4522")
plot(cors[[3]], main="df_4358")
plot(cors[[4]], main="df_4231")
plot(cors[[5]], main="df_4218")
plot(cors[[6]], main="df_4136")
plot(cors[[7]], main="df_3502")
plot(cors[[8]], main="df_2821")
plot(cors[[9]], main="df_1917")

i = 1
for (cor in cors) {
  cor_means[[i]] = mean(cors[[i]])
  i = i +1 
}

# > cor_means
################################################
# Notes:
# [1] 0.9662468
# [1] 0.9519713
# [1] 0.918355
# [1] 0.007753404
# [1] 0.9957909
# [1] 0.9364576
# [1] 0.053371
# [1] 0.939854
# [1] 0.9841476
# Let's drop 4 and 7 (no correlation betweeen samples - investigate?)





# Compute correlations within disease samples from each dataset
cors = list()
i = 1
for (df in disease_dfs) {
  message(i, "-----------------------------------------------")
  print(dim(df))
  str(df)
  # Transpose to use cor()
  df_t = t(df)
  # Compute corr
  corre = cor(df_t, use="complete.obs")
  cors[[i]] = corre
  i = i + 1
}
names(cors) = c('df_4523','df_4522','df_4358','df_4218','df_4136','df_2821','df_1917')
attach(cors)
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
plot(cors[[1]], main="df_4523")
plot(cors[[2]], main="df_4522")
plot(cors[[3]], main="df_4358")
plot(cors[[4]], main="df_4218")
plot(cors[[5]], main="df_4136")
plot(cors[[6]], main="df_2821")
plot(cors[[7]], main="df_1917")

i = 1
cor_means = list()
for (cor in cors) {
  cor_means[[i]] = mean(cors[[i]])
  i = i +1 
}

# > cor_means
################################################
# Notes:
# [1] 0.9682622
# [1] 0.9601459
# [1] 0.912433
# [1] 0.8739121
# [1] 0.9695674
# [1] 0.9458787
# [1] 0.9639092
# Good correlations among these dz samples



# Average across rows for each control state
avgs_ctrl = list()
i = 1
for (df in control_dfs) {
  message(i, " is [", nrow(df), " x ", ncol(df), "] Matrix --------------------")
  colMeanv = colMeans(df, na.rm=TRUE)
  # Having the same-name function and variable is not a good idea
  avgs_ctrl[[i]] = colMeanv
  i = i + 1
}

# Average across rows for each disease state
avgs_dz = list()
i = 1
for (df in disease_dfs) {
  message(i, " is [", nrow(df), " x ", ncol(df), "] Matrix --------------------")
  colMeanv = colMeans(df, na.rm=TRUE)
  avgs_dz[[i]] = colMeanv
  i = i + 1
}

# Subtract the control average expression values from the dz averages
normalized = list()
i = 1
print(length(avgs_ctrl))
print(length(avgs_dz))
for (i in 1:length(avgs_ctrl)) {
  message(i, " -----------------------------------------------")
  normalized[[i]] = avgs_dz[[i]] - avgs_ctrl[[i]]
  i = i + 1
}


require(reshape2)
# worked = list(4,5,7,8,9,10,11,12,13,14,15,16,17) # (also 2 worked)
# worked_2 = list(9,11,14,15,16)
# worked_3 = list(4,5,7,10,12,13,17)
# worked_4 = list(7,10)
# avg_dz_exp_vals = as.data.frame(avgs[7])
# avg_dz_exp_vals = data.frame(as.data.frame(avgs[7]))
# for (i in worked_4) {
#   print(i)
#   avg_dz_exp_vals <- cbind(avg_dz_exp_vals, as.data.frame(avgs[i]))
# }



# Averages to DF
normalized = as.data.frame(normalized)


###########################################################
# Compute correlation between diseases expression values  #
###########################################################
install.packages("corrplot")
library(corrplot)

list(df_4523,df_4522,df_4358,df_4218,df_4136,df_2821,df_1917)

<<<<<<< HEAD

=======
>>>>>>> ae9e79fcd6b937d94da196434cc28ab466d0baa5
# Again, R already has functions c() and names().
cv = cor(normalized, use="complete.obs") # Cor vector
namev = c("GDS4523 - \nSchizophrenia",
          "GDS4522 - \nSchizophrenia",
          "GDS4358 - \nHIV",
          "GDS4218 - \nMultiple Sclerosis",
          "GDS4136 - \nAlzheimers",
          "GDS2821 - \nParkinsons",
          "GDS1917 - \nSchizophrenia")
colnames(cv) <- rownames(cv) <- namev

<<<<<<< HEAD
attach(cv)

=======
attach(as.data.frame(c))
attach(cv)
>>>>>>> ae9e79fcd6b937d94da196434cc28ab466d0baa5
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
#corrplot.mixed(cv, t1.pos="r", t1.col="blue", c1.srt=60, c1.pos="r", cl.align.text="r", mar=c(1,1,1,1), height=1600, width=1600)
# perhaps you'd like to do this?
corrplot.mixed(cv, tl.pos=c("d", "lt", "n")[1], tl.col="blue", cl.align.text="r", mar=c(1,1,1,1))

# More highly correlatd than last time?? √

pData(ESET_3502)
exprs(ESET_3502)
assayDataElement(ESET_3502)

plot_scatterplot_matrix <- function(normalized){
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y, use = 'complete.obs')) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
    
    test <- cor.test(x,y) 
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")) 
    
    text(0.5, 0.5, txt, cex = cex *.8) 
    text(.8, .8, Signif, cex=cex, col=2) 
  }
  
  panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
                          cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    # if (any(ok)) 
    #     lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
    #           col = col.smooth, ...)
  }
  
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
  }
  
  pairs(normalized, lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist )
  
}
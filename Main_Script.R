all_GDSv <- c(5204,4879,4859,4838,4758,4532,4522,4523,4477,4414,4358,4231,4218,4154,4136,4135,3834,3502,3459,3345,3129,3128,3113,3110,3069,2978,2941,2821,2795,2613,2191,2190,2154,1962,1917,1912,1835,1816,1815,1813,1726,1253,1096,1085,910,909,833,810,707,596,564,426,232,181)

pfc = c(2190, 3502, 4414, 4523) #, 4532) - not in same region



# Download all the datasets
all_GDS = list()
j = 1
for (i in all_GDSv){
  all_GDS[[j]] = assign(paste('GDS', all_GDSv[i], sep='_'), getGEO(GEO=paste("GDS", i, sep='')))
  j = j + 1
}

# Extract all to ESETs
all_ESET = list()
j = 1
for (gds in all_GDS){
  all_ESET[[j]] = assign(paste('ESET', all_GDSv[j], sep='_'), GDS2eSet(GDS=gds, do.log2 = TRUE)) 
  j = j + 1
}

# To df's
all_DF = list()
j = 1;
for (eset in all_ESET) { 
  all_DF[[j]] = assign(paste('df', all_GDSv[j], sep='_'), as.data.frame(eset))
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

# Find NAs, NANs (NANs are a todo)
nan_count <-function (x) sapply(x, function(y) sum(is.nan(y)))
for (df in gpl570) {
  print(any(is.na(df)))
  # print(any(is.nan(df)))
}


#####################################################################
# Notes:

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
for (df in gpl570) {
  common_cols <- intersect(common_cols, colnames(df))
}

# --> 54677 Common Columns!


# Subset datasets to just disease states

# Compute correlations within samples from each dataset
disease_dfs = list()
i = 1
for (df in gpl570) {
  print(i)
  print("-----------------------------------------------")
  str(df)
  # Subset to just DZ instances
  df_disease = df[,df$disease.state != "control"]
  # Remove uncommon columns
  df_disease = df_disease[, (names(df_disease) %in% common_cols)]
  # Remove metadata
  drops <- c("sample", "age", "gender", "tissue", "genotype.variation", "development.stage", "agent", "other", "cell.type", "disease.state", "description")
  df_disease = df_disease[, !(names(df_disease) %in% drops)]
  disease_dfs[[i]] = df_disease
  i = i +1
}



# Compute correlations within samples from each dataset
cors = list()
i = 1
for (df in disease_dfs) {
  print(i)
  print("-----------------------------------------------")
  str(df)
  # Transpose to use cor()
  df_disease_t = t(df)
  # Compute corr
  corre = cor(df_disease_t)
  cors[[i]] = corre
  i = i + 1
}

names(cors) = c('df_5204','df_4838','df_4532','df_4523','df_4522','df_4477','df_4358','df_4231','df_4218','df_4136','df_4135','df_3502','df_2821','df_2795','df_2154','df_1962','df_1917')
attach(cors)
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
plot(cors[[2]], main="asdf")
plot(cors[[4]], main="asdf")
plot(cors[[7]], main="asdf")
plot(cors[[10]], main="asdf")
plot(cors[[14]], main="asdf")
plot(cors[[15]], main="asdf")
plot(cors[[17]], main="asdf")



# Average across rows for each disease state
avgs = list()
i = 1
for (df in disease_dfs) {
  print(i)
  print(dim(df))
  print("-----------------------------------------------")
  colMeans = colMeans(df, na.rm=TRUE)
  avgs[[i]] = colMeans
  i = i + 1
}



require(reshape2)
worked = list(4,5,7,8,9,10,11,12,13,14,15,16,17) # (also 2 worked)
worked_2 = list(9,11,14,15,16)
avg_dz_exp_vals = as.data.frame(avgs[2])
for (i in worked_2) {
  avg_dz_exp_vals <- cbind(avg_dz_exp_vals, as.data.frame(avgs[i]))
}



install.packages("corrplot")
library(corrplot)
# Compute correlation between diseases expression values
c = cor(avg_dz_exp_vals, use="complete.obs")
names = c("GDS4838 - Neuroectodermal Tumors",
          "GDS4218 - Multiple Sclerosis Brain Lesions",
          "GDS4135 - Temporal Cortex, Ageing Study",
          "GDS2795 - Alzheimer's Disease",
          "GDS2154 - Inflammatory Dilated Cardiomyopathy",
          "GDS1962 - Glioma tissue")
colnames(c) <- rownames(c) <- names
corrplot.mixed(c, t1.pos="r", t1.col="blue", c1.srt=60, c1.pos="r", cl.align.text="r", mar=c(1,1,1,1), height=1600, width=1600)




pData(ESET_3502)
exprs(ESET_3502)
assayDataElement(ESET_3502)


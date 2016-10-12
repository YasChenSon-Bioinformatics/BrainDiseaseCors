#################################################################################################################################
# GDS1917, GPL570, Cerebellar cortex, Cerebellar cortex in schizophrenia

# Subset to just DZ instances
df_1917_disease = df_1917[,df_1917$disease.state != "control"]

# Remove metadata
drops <- c("sample", "disease.state", "description")
df_1917_disease = df_1917_disease[ , !(names(df_1917_disease) %in% drops)]

# Transpose to use cor()
df_1917_disease_t = t(df_1917_disease)

# Compute corr
cor_1917 = cor(df_1917_disease_t)



#################################################################################################################################
# GDS4136, GPL570, hippocampal CA1 gray matter, Various stages of Alzheimer's disease: laser-captured hippocampal CA1 gray matter

# Subset to just DZ instances
df_4136_disease = df_4136[,df_4136$disease.state != "control"] # Other values are "incipient stage", "moderate stage", "severe stage"

# Remove metadata
drops <- c("sample", "age", "disease.state", "description")
df_4136_disease = df_4136_disease[ , !(names(df_4136_disease) %in% drops)]

# Transpose to use cor()
df_4136_disease_t = t(df_4136_disease)

# Compute corr
cor_4136 = cor(df_4136_disease_t)

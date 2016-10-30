
library(GEOquery) # installed by BiocInstaller::biocLite("GEOQuery")
library(limma)

   GDS4523 <- getGEO('GDS4523')

  ESET4523 <- GDS2eSet(GDS4523, do.log2 = TRUE)  # The original GDS data is not log transformed.

MATRIX4523 <- as.matrix( ESET4523 ) # call as.matrix.ExpressionSet defined in library(limma)

index_control <- GDS4523@dataTable@columns$disease.state == 'control'
index_disease <- GDS4523@dataTable@columns$disease.state != 'control'
reported_probes <- c('200011_s_at') # 1 out of 49 in http://www.nature.com/mp/journal/v14/n12/fig_tab/mp200918t2.html#figure-title

gene_expression_values <- exprs(ESET4523) # extract gene expression values

# Analysis 1. two sample t-test
# (Probe '200011_s_at' expression value difference between the mean of 'control' samples and the mean of 'schizophrenia' samples)

ttested <- 
t.test( gene_expression_values[reported_probes, index_control], var.equal = FALSE,
        gene_expression_values[reported_probes, index_disease], alternative = 'two.sided')
# p-value 0.256 (t=1.1504). Not the same as the reported p-values.

# Analysis 2. linear model with empirical bayes

lmfitted <- lmFit(MATRIX4523, design = model.matrix( ~ GDS4523@dataTable@columns$disease.state))
ebayesd <- eBayes(lmfitted)                 # not ebayes (backward-compatibility) but eBayes
topped <- topTable(ebayesd, number = 60000) # not toptable (backward-compatibility) but topTable

topped[rownames(topped) == reported_probes, ]
# Unadjustd P.Value is 0.1960. Adjusted 0.9999. Not the same as the reported p-values.

# Analysis 3. Diagnosis by visual plots

mean_control <- mean(gene_expression_values[reported_probes, index_control])
mean_disease <- mean(gene_expression_values[reported_probes, index_disease])

plot(gene_expression_values[reported_probes, ], col=as.factor(index_control), pch=as.numeric(factor(gender)), type='n',
     main = paste0("Gene Expression Values (log2-transformed) of probe ", reported_probes," in GDS4523 ",
                   "\n Black : Schizophrenia samples (28)    Red: Control samples (23) \n",
                   " p-value of 2 sample t-test :", round(ttested$p.value,4) ))
text(gene_expression_values[reported_probes, ], labels=age)
segments(x0=min(which(index_control==TRUE)) , y0=mean_control,
         x1=max(which(index_control==TRUE)) , y1=mean_control, col="red")

segments(x0=min(which(index_disease==TRUE)) , y0=mean_disease,
         x1=max(which(index_disease==TRUE)) , y1=mean_disease, col="black")

# Not so different


# Analysis 4. Consideration of gender, age

dz <- GDS4523@dataTable@columns$disease.state
gender <- GDS4523@dataTable@columns$gender
age <- as.numeric(substr(GDS4523@dataTable@columns$age,1,2))
lmfitted2 <- lmFit(MATRIX4523, design = model.matrix( ~  dz + gender + age))
ebayesd2 <- eBayes(lmfitted2)                 # not ebayes (backward-compatibility) but eBayes
topped2 <- topTable(ebayesd2, number = 60) # not toptable (backward-compatibility) but topTable

topped2[rownames(topped2) == reported_probes, ]

plot(topped2)

plot(age, gene_expression_values[reported_probes, ])
hist(gene_expression_values[reported_probes, gender=='male'])
hist(gene_expression_values[reported_probes, gender!='male'])
text(gene_expression_values[reported_probes, ], labels=age, col=as.numeric(gender))

# Analysis 5.
i <- 2
dummy <- '45288_at'
material <- exprs(ESETl[[i]])[ rownames(ESETl[[i]]) %in% c(reported_probes, dummy), ] %>% as.data.frame(.,stringsAsFactors=FALSE) %>% rownames_to_column('probe') %>% gather(sample,eval, -starts_with('probe')) %>% mutate( sample = factor(sample)) %>% left_join(., GDSl[[i]]@dataTable@columns, by='sample') %>% mutate( age = gsub('[^0-9]','',age) %>% as.numeric, meta = paste0(age,'',substr(gender,1,1)) )

source('~/BrainDiseaseCors/caprice/R/pairs-extension.R')
library(GGally)
m <- material %>% filter( probe == reported_probes) %>% dplyr::select(age,gender,disease.state,eval) %>% mutate( disease.state=factor(disease.state) %>% as.numeric%>% jitter(), gender = as.factor(gender) %>% as.numeric %>% jitter())
#pairs(, lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
ggpairs(m, title = paste0("GDS4523 SHOULD BE ANALYZED WITH AGE CONSIDERATION\n ", reported_probes, " (reported in paper)"), diag = list(discrete="barDiag", continuous='barDiag'))


ggplot(material,aes(x=sample,y=eval,color=disease.state)) +  geom_text(aes(label=meta),size=3,lineheight=.6) + facet_wrap(  ~  strip, ncol=5) + theme(legend.position="bottom", axis.text.x = element_text(angle=90)) + labs(title = paste0(GDSl[[i]]@header$title,'  -  ', GDSl[[i]]@header$dataset_id[1]), subtitle = strwrap(GDSl[[i]]@header$description[1], width = 140) %>% paste(.,collapse='\n') )

a <- 
  exprs(ESETl[[i]])[ rownames(ESETl[[1]]) %in% c('45288_at','38703_at'), ] %>% as.data.frame(.,stringsAsFactors=FALSE) %>% rownames_to_column('probe') %>% gather(sample,eval, -starts_with('probe')) %>% mutate( sample = factor(sample)) %>% left_join(., GDSl[[i]]@dataTable@columns, by='sample')

t.test( (a %>% filter(disease.state=='control' & probe != '38703_at'))$eval,
        (a %>% filter(disease.state!='control' & probe != '38703_at' ))$eval )

ggplot(a,aes(x=sample,y=eval,color=disease.state)) +  geom_text(aes(label=probe),size=3,lineheight=.6)+ facet_wrap(  ~  probe, ncol=5)


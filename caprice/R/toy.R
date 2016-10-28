

dogdf <- data.frame( sibship = c(1,1,2,2,3,3),
                     treat   = c('C','T','C','T','C','T'))


design <- model.matrix( ~ factor(dogdf$sibship) + factor(dogdf$treat) )

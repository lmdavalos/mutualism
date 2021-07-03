##libraries
library(brms)

##remove prior data
rm(list=ls()) 

##get data
## + 1
load("multiresponse_v3.RData")

##first make list with relevant models
mylist<-list(m.car.b, m.cym.p, m.ger.d, m.hep.2, m.myr.b, m.non.2, m.phe.a, m.phe.b)

##rename every model to be compared and add loo criterion at the same time
p_one<-sapply(mylist, function(mylist) add_criterion(mylist, "loo"), simplify=F)

## + 0.1
load("multiresponse_v3_0.RData")

p_pon<-sapply(mylist, function(mylist) add_criterion(mylist, "loo"), simplify=F)

## + 10
load("multiresponse_v3_10.RData")

p_ten<-sapply(mylist, function(mylist) add_criterion(mylist, "loo"), simplify=F)

##compare loo for each model
coml<-list()

for(i in 1:length(mylist))
{	coml[[i]]<-loo_compare(p_one[[i]], p_pon[[i]], p_ten[[i]], criterion = "loo")}

##print out comparisons
##all models are EXACTLYthe same

sink("model_comparison.txt")
sapply(coml, function(coml) print(coml))
sink()

##save
save.image("compare.RData")

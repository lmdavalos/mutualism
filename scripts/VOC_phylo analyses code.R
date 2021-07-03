library(parallelsugar)
library(geiger)
library(phytools)
library(mvMORPH)
library(l1ou)


mds.scores<-read.csv(file="Piper_MDSscores.csv",row.names=1)
raw.VOC<-read.csv(file="ripefruit_VOCaver.csv",row.names=1)
total.VOC<-apply(raw.VOC,1,sum)


###### time calibrating the phylogeny with chronos
###### tree may be slightly different, set.seed not set, tree file provided

VOCtree<-read.nexus(file="VOCtree.nex")
VOCtree<-root(VOCtree,resolve.root=TRUE,outgroup=VOCtree$tip.label[17:22])
VOCtree$edge.length[1]<-0.00000001

cals<-makeChronosCalib(VOCtree)
cals[1,]<-c(36,54.4,56.1,FALSE)
cals[2,]<-c(23,67,72,FALSE) #neotropical clade maximum

time.tree<-chronos(VOCtree,calibration=cals)
class(time.tree)<-"phylo"


### plotting traits on tree
par(mfrow=c(6,2))
for(i in 1:6){
x<-mds.scores[,i];names(x)<-row.names(mds.scores)
contMap(time.tree,x)
phenogram(time.tree,x)
}

physignal(mds.scores[,1:3],time.tree)

sig.mat<-matrix(ncol=2,nrow=6)
for(i in 1:6){
  x<-mds.scores[,i];names(x)<-row.names(mds.scores)
 sig<-phylosig(time.tree,x,test=TRUE)
 sig.mat[i,]<-c(sig$K,sig$P)
 }



##### model fitting multivariate BM/OU/EB models for mds.scores
BM<-mvBM(time.tree,mds.scores[,1:3])
OU<-mvOU(time.tree,mds.scores[,1:3])
EB<-mvEB(time.tree,mds.scores[,1:3])


############### running l1ou on the mds.scores
ltree<-adjust_data(time.tree,mds.scores)$tree
VOC.test<-estimate_shift_configuration(ltree,mds.scores[ltree$tip.label,1:3],criterion="pBIC",nCores=4,)


############ running l1ou on simulated data (BM)
  rm<-ratematrix(time.tree,total.VOC)
  sims<-sim.char(time.tree,rm,model="BM",nsim=100)
  
  res.list<-vector(length=100,mode="list")
  
  for(i in 1:100){
    progress(i,max.value=100)
    sim.test<-NA
    sim.cons<-NA
    sim.test<-try(estimate_shift_configuration(ltree,sims[ltree$tip.label,1,i],criterion="pBIC",nCores=4,rescale=FALSE),silent=TRUE)
    res.list[[i]]<-sim.test
  }  



######## cophylogenetic comparison of piper phylogeny and VOCs
d<-vegdist(raw.VOC,"bray")
fit<-hclust(d,method="average")
pfit<-as.phylo(fit)
cops<-cophylo(time.tree,pfit)



######### preferred VOCs as discrete
voc.VIP<-read.csv(file="ripefruit_aver_prefVOCs.csv")
vocd<-ifelse(voc.VIP[,3:6]>0,1,0)
row.names(vocd)<-voc.VIP[,1]


s.car<-make.simmap(time.tree,vocd[,1],model="ARD",nsim=100)
s.x2hep<-make.simmap(time.tree,vocd[,2],model="ARD",nsim=100) 
s.germ<-make.simmap(time.tree,vocd[,3],model="ARD",nsim=100) 
s.bfar<-make.simmap(time.tree,vocd[,4],model="ARD",nsim=100) 




### fit models for transition rates
fit.ARD.car<-fitMk(time.tree,vocd[,1],model="ARD")
fit.ARD.x2hep<-fitMk(time.tree,vocd[,2],model="ARD")
fit.ARD.germ<-fitMk(time.tree,vocd[,3],model="ARD")
fit.ARD.bfar<-fitMk(time.tree,vocd[,4],model="ARD")

fit.ER.car<-fitMk(time.tree,vocd[,1],model="SYM")
fit.ER.x2hep<-fitMk(time.tree,vocd[,2],model="ER")
fit.ER.germ<-fitMk(time.tree,vocd[,3],model="SYM")
fit.ER.bfar<-fitMk(time.tree,vocd[,4],model="SYM")




#running l1ou on log transformed VOC data
vocs<-voc.VIP[,3:6];row.names(vocs)<-voc.VIP[,1]
vocs2<-log10(vocs+1)
ltree<-adjust_data(time.tree,voc2)$tree
voc.test<-estimate_shift_configuration(ltree,vocs2[time.tree$tip.label,2],criterion="pBIC",rescale=FALSE)



###simulations for l1ou
rm<-ratematrix(time.tree,vocs2)
sims<-sim.char(time.tree,rm,model="BM",nsim=100)

res.list<-vector(length=100,mode="list")

for(i in 1:100){
  progress(i,max.value=100)
  sim.test<-NA
  sim.cons<-NA
  sim.test<-try(estimate_shift_configuration(ltree,sims[ltree$tip.label,,i],criterion="pBIC",rescale=FALSE),silent=TRUE)
  res.list[[i]]<-sim.test
}  

aces<-apply(vocs2,2,function(x){ace(x,time.tree)$ace[1]})


sims<-array(dim=c(22,5,100))
for(i in 1:5){
  sims[,i,]<-apply(sims,3,function(x){fastBM(tree=time.tree,sig2=rm[i,i],bounds=c(0,Inf),a=aces[i])})
}
row.names(sims)<-time.tree$tip.label






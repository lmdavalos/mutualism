library(ape)

##remove prior data
rm(list=ls()) 

##get the data

diet<-read.csv("Carollia_pruneddiet_correct_names.csv")
mds<-read.csv("Piper_MDSscores.csv")
voc<-read.csv("Piper_ripefruit_VOCmeans_total_shannon_abundance.csv")
##ind<-read.csv("VOC_ripefruit_perindiv.csv")

##get  the tree
tree<-read.tree("time.tree.correct.names2.phy")

##get phylogenetic matrix
A <- ape::vcv.phylo(tree)

##get the data together
dat1<-merge(diet,mds)

##proportions from 0-1
dat1$cas <-dat1$Carollia_castanea/100
dat1$per <-dat1$Carollia_perspicillata/100
dat1$sow <-dat1$Carollia_sowelli/100

#individual variation cannot be used bc of lack of matching observations
##get the data together averages
dat2<-merge(diet,voc)

##proportions from 0-1
dat2$cas <-dat2$Carollia_castanea/100
dat2$per <-dat2$Carollia_perspicillata/100
dat2$sow <-dat2	$Carollia_sowelli/100

##prepare abundance data
tr1<-drop.tip(tree, c("Piper_aduncum", "Piper_peltatum", "Piper_umbricola"))

##get new phylogenetic matrix
As <- ape::vcv.phylo(tr1)

save.image("setup_VOCs.RData")

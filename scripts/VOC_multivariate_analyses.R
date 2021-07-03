# Multivariate analyses of Piper VOCs
library(vegan)
library(fpc)
library(cluster)

#Import data
df<-read.csv("Piper_ripefruit_allVOCmeans.csv", row.names = 'Species')
row.names(df)
y<-df
dim(y)

#MDS analysis
mds<-capscale(y~1, distance = "bray")
summary(mds)

sink("Results_VOC_MDS.txt")
print(summary(mds))
sink()

mds.scores<-scores(mds, choices=c(1:21), display=c("wa"), scaling=2)
write.csv(mds.scores, file = "Piper_MDSscores.csv")

#Plot
plot(mds,display='sites',type='none', xlim = c(-2,2), ylim = c(-2,2),choices=c(1,2))
points(mds,choices=c(3,4),'sites',pch=19 ,cex=1)
text(mds, labels = row.names(df),cex = 1,pos=2, choices=c(3,4))

#add vectors
ef15<-envfit(pcoa,y,permu=1000,choices=c(1:5))
ef15$vectors
plot(ef15,col='blue', p.max=.001, cex = 0.8, choices=c(3,4))

options(max.print=999999)

sink("Vectors_VOC_MDS.txt")
print(ef15)
sink()

#calculating Shannon diversity index for VOCs
div <- diversity(df, index="shannon")
write.csv(div, file = "VOC_shannon.csv")

##libraries
library(brms)
library(ggplot2)

##remove prior data
rm(list=ls()) 

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##set theme ggplot2
theme_set(theme_minimal())

##get data
## + 1
load("multiresponse_v3.RData")

##first make list with relevant models
mylist<-list(m.abu.s, m.shw.e, m.tot.e, m.num.e, 
m.car.b, m.cym.p, m.ger.d, m.hep.2, m.myr.b, m.non.2, m.phe.a, m.phe.b, 
m.mds.1, m.mds.2, m.mds.3, m.mds.4, m.mds.5,
m.hep.2, m.hep.2, m.hep.2)

##make vector with corresponding variable
myvec<-c("abundancesalazar", "shannon", "log10totalemission", "numvocs",
"log10Caryophyllene1", "log10pCymene1", "log10GermacreneD1", "log10X2Heptanol1", "log10betaMyrcene1", "log10X2Nonanone1", "log10alphaPhellandrene1", "log10betaPhellandrene1",
"MDS1", "MDS2", "MDS3", "MDS4", "MDS5",
"cas", "per", "sow")

#generate pic for each model loo-pit
pplo<-list()

for(i in 1:length(mylist))
{	pplo[[i]]<-pp_check(mylist[[i]], resp=myvec[[i]], type = "loo_pit_overlay" )}

##print loopit
pdf("voc_stat_loop.pdf", h=3, w=12)
multiplot(pplo[[4]], pplo[[2]], pplo[[3]], pplo[[1]], cols=4)
dev.off()

##print
pdf("vocs_indi_loop.pdf", h=6, w=15)
multiplot(pplo[[8]], pplo[[9]], pplo[[10]], pplo[[12]], pplo[[11]], pplo[[7]], pplo[[5]], pplo[[6]], cols=4)
dev.off()

##print
pdf("mds_loop.pdf", h=3, w=15)
multiplot(pplo[[13]], pplo[[14]], pplo[[15]], pplo[[16]], pplo[[17]], cols=5)
dev.off()

##print
pdf("spp_loop.pdf", h=3, w=15)
multiplot(pplo[[18]], pplo[[19]], pplo[[20]], cols=3)
dev.off()

#generate pic for each model ppcheck loo ribbon
ppl1<-list()

for(i in 1:length(mylist))
{	ppl1[[i]]<-pp_check(mylist[[i]], resp=myvec[[i]], type = "loo_ribbon" )}

##print loopit
pdf("voc_stat_ppch.pdf", h=3, w=12)
multiplot(ppl1[[4]], ppl1[[2]], ppl1[[3]], ppl1[[1]], cols=4)
dev.off()

##print
pdf("vocs_indi_ppch.pdf", h=6, w=15)
multiplot(ppl1[[8]], ppl1[[9]], ppl1[[10]], ppl1[[12]], ppl1[[11]], ppl1[[7]], ppl1[[5]], ppl1[[6]], cols=4)
dev.off()

##print
pdf("mds_ppch.pdf", h=3, w=15)
multiplot(ppl1[[13]], ppl1[[14]], ppl1[[15]], ppl1[[16]], ppl1[[17]], cols=5)
dev.off()

##print
pdf("spp_ppch.pdf", h=3, w=15)
multiplot(ppl1[[18]], ppl1[[19]], ppl1[[20]], cols=3)
dev.off()

#generate pic for each model ppcheck ribbon
ppl2<-list()

for(i in 1:length(mylist))
{	ppl2[[i]]<-pp_check(mylist[[i]], resp=myvec[[i]], type = "ribbon" )}

##print ribbon
pdf("voc_stat_ribb.pdf", h=3, w=12)
multiplot(ppl2[[4]], ppl2[[2]], ppl2[[3]], ppl2[[1]], cols=4)
dev.off()

##print
pdf("vocs_indi_ribb.pdf", h=6, w=15)
multiplot(ppl2[[8]], ppl2[[9]], ppl2[[10]], ppl2[[12]], ppl2[[11]], ppl2[[7]], ppl2[[5]], ppl2[[6]], cols=4)
dev.off()

##print
pdf("mds_ribb.pdf", h=3, w=15)
multiplot(ppl2[[13]], ppl2[[14]], ppl2[[15]], ppl2[[16]], ppl2[[17]], cols=5)
dev.off()

##print
pdf("spp_ribb.pdf", h=3, w=15)
multiplot(ppl2[[18]], ppl2[[19]], ppl2[[20]], cols=3)
dev.off()

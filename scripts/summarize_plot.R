##library
library(brms)
library(ggplot2)
library(bayesplot)

##clear out cache
rm(list=ls())

##load data
load("multiresponse_v3.RData")

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


##print summary MDS
sink("mds_coef.txt")
print("MDS 1")
print(summary(m.mds.1))

print("MDS 2")
print(summary(m.mds.2))

print("MDS 3")
print(summary(m.mds.3))

print("MDS 4")
print(summary(m.mds.4))

print("MDS 5")
print(summary(m.mds.5))
sink()

##print summary emissions
##print summary
sink("vocs_stat_coef.txt")
print("number of VOCs")
print(summary(m.num.e))

print("tota emissions")
print(summary(m.tot.e))

print("Shannon Wiener")
print(summary(m.shw.e))

print("Piper abundance")
print(summary(m.abu.s))
sink()

##print individual VOCs
sink("vocs_indi_coef.txt")
print("X2.Heptanol")
print(summary(m.hep.2))

print("X2.Nonanone")
print(summary(m.non.2))

print("alpha.Phellandrene")
print(summary(m.phe.a))

print("beta.Caryophyllene")
print(summary(m.car.b))

print("beta.Myrcene")
print(summary(m.myr.b))

print("beta.Phellandrene")
print(summary(m.phe.b))

print("p.Cymene")
print(summary(m.cym.p))

print("Germacrene.d")
print(summary(m.ger.d))
sink()

##make less terrible
color_scheme_set("gray")
bayesplot_theme_set(theme_minimal())

##MDS plots
firs<-mcmc_areas(m.mds.1, prob=.9, area_method="equal height", pars=c("b_cas_MDS1", "b_per_MDS1", "b_sow_MDS1"), prob_outer=0.95) +
 ggplot2::labs(x = "MDS1")

seco<-mcmc_areas(m.mds.2, prob=.9, area_method="equal height", pars=c("b_cas_MDS2", "b_per_MDS2", "b_sow_MDS2"), prob_outer=0.95) +
 ggplot2::labs(x = "MDS2")

thir<-mcmc_areas(m.mds.3, prob=.9, area_method="equal height", pars=c("b_cas_MDS3", "b_per_MDS3", "b_sow_MDS3"), prob_outer=0.95) +
 ggplot2::labs(x = "MDS3")

four<-mcmc_areas(m.mds.4, prob=.9, area_method="equal height", pars=c("b_cas_MDS4", "b_per_MDS4", "b_sow_MDS4"), prob_outer=0.95) +
 ggplot2::labs(x = "MDS4")

five<-mcmc_areas(m.mds.5, prob=.9, area_method="equal height", pars=c("b_cas_MDS5", "b_per_MDS5", "b_sow_MDS5"), prob_outer=0.95) +
 ggplot2::labs(x = "MDS5")

##print
pdf("mds_coefs.pdf", h=3, w=15)
multiplot(firs, seco, thir, four, five, cols=5)
dev.off()

##VOC stats and abundance plots
numv<-mcmc_areas(m.num.e, prob=.9, area_method="equal height", regex_pars="_num_vocs", prob_outer=0.95) +
 ggplot2::labs(x = "Count of VOCs")

tote<-mcmc_areas(m.tot.e, prob=.9, area_method="equal height", regex_pars="log10total_emission", prob_outer=0.95) +
 ggplot2::labs(x = "Total emissions (log10)")

shan<-mcmc_areas(m.shw.e, prob=.9, area_method="equal height", pars=c("b_cas_shannon", "b_per_shannon", "b_sow_shannon"), prob_outer=0.95) +
 ggplot2::labs(x = "Shannon-Wiener")

abun<-mcmc_areas(m.abu.s, prob=.9, area_method="equal height", pars=c("b_cas_abundance_salazar", "b_per_abundance_salazar", "b_sow_abundance_salazar"), prob_outer=0.95) +
 ggplot2::labs(x = "Abundance in field (count)")

 
##print
pdf("voc_stat_coefs.pdf", h=3, w=12)
multiplot(numv, shan, tote, abun, cols=4)
dev.off()

##individual VOCs
hep2<-mcmc_areas(m.hep.2, prob=.9, area_method="equal height", pars=c("b_cas_log10X2.HeptanolP1", "b_per_log10X2.HeptanolP1", "b_sow_log10X2.HeptanolP1"), prob_outer=0.95) +
 ggplot2::labs(x = "2-heptanol")

non2<-mcmc_areas(m.non.2, prob=.9, area_method="equal height", pars=c("b_cas_log10X2.NonanoneP1", "b_per_log10X2.NonanoneP1", "b_sow_log10X2.NonanoneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "2-nonanone")

aphe<-mcmc_areas(m.phe.a, prob=.9, area_method="equal height", pars=c("b_cas_log10alpha.PhellandreneP1", "b_per_log10alpha.PhellandreneP1", "b_sow_log10alpha.PhellandreneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "alpha-phellandrene")

bcar<-mcmc_areas(m.car.b, prob=.9, area_method="equal height", pars=c("b_cas_log10CaryophylleneP1", "b_per_log10CaryophylleneP1", "b_sow_log10CaryophylleneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "beta-caryophyllene")

bmyr<-mcmc_areas(m.myr.b, prob=.9, area_method="equal height", pars=c("b_cas_log10beta.MyrceneP1", "b_per_log10beta.MyrceneP1", "b_sow_log10beta.MyrceneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "beta-myrcene")

bphe<-mcmc_areas(m.phe.b, prob=.9, area_method="equal height", pars=c("b_cas_log10beta.PhellandreneP1", "b_per_log10beta.PhellandreneP1", "b_sow_log10beta.PhellandreneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "beta-phellandrene")

gerd<-mcmc_areas(m.ger.d, prob=.9, area_method="equal height", pars=c("b_cas_log10Germacrene.DP1", "b_per_log10Germacrene.DP1", "b_sow_log10Germacrene.DP1"), prob_outer=0.95) +
 ggplot2::labs(x = "germacrene-d")

pcym<-mcmc_areas(m.cym.p, prob=.9, area_method="equal height", pars=c("b_cas_log10p.CymeneP1", "b_per_log10p.CymeneP1", "b_sow_log10p.CymeneP1"), prob_outer=0.95) +
 ggplot2::labs(x = "p-cymene")
 
##print
pdf("vocs_indi_coefs.pdf", h=6, w=15)
multiplot(hep2, bmyr, non2, bphe, aphe, gerd, bcar, pcym, cols=4)
dev.off()


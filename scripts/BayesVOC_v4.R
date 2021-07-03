##library
library(brms)

##clear out cache
rm(list=ls())

##load data
load("setup_VOCs.RData")

##you may need to increase iterations, raise adapt_delta and treedepth and rerun to eliminate or reduce divergent transitions 

##X2.Nonanone
non.2 <- bf(log10(X2.Nonanone+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.o <- bf(cas ~ log10(X2.Nonanone+1) + (1|gr(Species, cov = A)))
per.o <- bf(per ~ log10(X2.Nonanone+1) + (1|gr(Species, cov = A)))
sow.o <- bf(sow ~ log10(X2.Nonanone+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.non.2 <- brm(non.2 +
        	cas.o +
        	per.o +
        	sow.o +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##beta.Myrcene
myr.b <- bf(log10(beta.Myrcene+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.m <- bf(cas ~ log10(beta.Myrcene+1) + (1|gr(Species, cov = A)))
per.m <- bf(per ~ log10(beta.Myrcene+1) + (1|gr(Species, cov = A)))
sow.m <- bf(sow ~ log10(beta.Myrcene+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.myr.b <- brm(myr.b +
        	cas.m +
        	per.m +
        	sow.m +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##X2.Heptanol
hep.2 <- bf(log10(X2.Heptanol+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.h <- bf(cas ~ log10(X2.Heptanol+1) + (1|gr(Species, cov = A)))
per.h <- bf(per ~ log10(X2.Heptanol+1) + (1|gr(Species, cov = A)))
sow.h <- bf(sow ~ log10(X2.Heptanol+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.hep.2 <- brm(hep.2 +
        	cas.h +
        	per.h +
        	sow.h +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##beta.Phellandrene
phe.b <- bf(log10(beta.Phellandrene+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.b <- bf(cas ~ log10(beta.Phellandrene+1) + (1|gr(Species, cov = A)))
per.b <- bf(per ~ log10(beta.Phellandrene+1) + (1|gr(Species, cov = A)))
sow.b <- bf(sow ~ log10(beta.Phellandrene+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.phe.b <- brm(phe.b +
        	cas.b +
        	per.b +
        	sow.b +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##alpha.Phellandrene
phe.a <- bf(log10(alpha.Phellandrene+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.a <- bf(cas ~ log10(alpha.Phellandrene+1) + (1|gr(Species, cov = A)))
per.a <- bf(per ~ log10(alpha.Phellandrene+1) + (1|gr(Species, cov = A)))
sow.a <- bf(sow ~ log10(alpha.Phellandrene+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.phe.a <- brm(phe.a +
        	cas.a +
        	per.a +
        	sow.a +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##S-W emissions is Gaussian
shw.e <- bf(shannon ~ (1|gr(Species, cov = A)))

##species representation
cas.e <- bf(cas ~ shannon + (1|gr(Species, cov = A)))
per.e <- bf(per ~ shannon + (1|gr(Species, cov = A)))
sow.e <- bf(sow ~ shannon + (1|gr(Species, cov = A)))

##run all models together at a time
m.shw.e <- brm(shw.e +
        	cas.e +
        	per.e +
        	sow.e +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##number emissions is Poisson, fancy!
num.e <- bf(num_vocs ~ (1|gr(Species, cov = A)))

##species representation
cas.n <- bf(cas ~ num_vocs + (1|gr(Species, cov = A)))
per.n <- bf(per ~ num_vocs + (1|gr(Species, cov = A)))
sow.n <- bf(sow ~ num_vocs + (1|gr(Species, cov = A)))

##run all models together at a time
m.num.e <- brm(num.e +
        	cas.n +
        	per.n +
        	sow.n +
            set_rescor(FALSE), 
            family = list("poisson",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 100000, warmup = 20000, thin = 200, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##VOC emissions
##total emission
tot.e <- bf(log10(total_emission) ~ (1|gr(Species, cov = A)))

##species representation
cas.t <- bf(cas ~ log10(total_emission) + (1|gr(Species, cov = A)))
per.t <- bf(per ~ log10(total_emission) + (1|gr(Species, cov = A)))
sow.t <- bf(sow ~ log10(total_emission) + (1|gr(Species, cov = A)))

##run all models together at a time
m.tot.e <- brm(tot.e +
        	cas.t +
        	per.t +
        	sow.t +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 10000, thin = 50, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##Germacrene.D
ger.d <- bf(log10(Germacrene.D+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.g <- bf(cas ~ log10(Germacrene.D+1) + (1|gr(Species, cov = A)))
per.g <- bf(per ~ log10(Germacrene.D+1) + (1|gr(Species, cov = A)))
sow.g <- bf(sow ~ log10(Germacrene.D+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.ger.d <- brm(ger.d +
        	cas.g +
        	per.g +
        	sow.g +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##model for MDS5
mds.5 <- bf(MDS5 ~ (1|gr(Species, cov = A)))

##species representation
cas.5 <- bf(cas ~ MDS5 + (1|gr(Species, cov = A)))
per.5 <- bf(per ~ MDS5 + (1|gr(Species, cov = A)))
sow.5 <- bf(sow ~ MDS5 + (1|gr(Species, cov = A)))

##run all models together at a time
m.mds.5 <- brm(mds.5 +
        	cas.5 +
        	per.5 +
        	sow.5 +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 5000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat1,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##model for MDS4
mds.4 <- bf(MDS4 ~ (1|gr(Species, cov = A)))

##species representation
cas.4 <- bf(cas ~ MDS4 + (1|gr(Species, cov = A)))
per.4 <- bf(per ~ MDS4 + (1|gr(Species, cov = A)))
sow.4 <- bf(sow ~ MDS4 + (1|gr(Species, cov = A)))

##run all models together at a time
m.mds.4 <- brm(mds.4 +
        	cas.4 +
        	per.4 +
        	sow.4 +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat1,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##model for MDS3
mds.3 <- bf(MDS3 ~ (1|gr(Species, cov = A)))

##species representation
cas.3 <- bf(cas ~ MDS3 + (1|gr(Species, cov = A)))
per.3 <- bf(per ~ MDS3 + (1|gr(Species, cov = A)))
sow.3 <- bf(sow ~ MDS3 + (1|gr(Species, cov = A)))

##run all models together at a time
m.mds.3 <- brm(mds.3 +
        	cas.3 +
        	per.3 +
        	sow.3 +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat1,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##model for MDS1
mds.1 <- bf(MDS1 ~ (1|gr(Species, cov = A)))

##species representation
cas.1 <- bf(cas ~ MDS1 + (1|gr(Species, cov = A)))
per.1 <- bf(per ~ MDS1 + (1|gr(Species, cov = A)))
sow.1 <- bf(sow ~ MDS1 + (1|gr(Species, cov = A)))

##run all models together at a time
m.mds.1 <- brm(mds.1 +
        	cas.1 +
        	per.1 +
        	sow.1 +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 5000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat1,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##p.Cymene
cym.p <- bf(log10(p.Cymene+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.p <- bf(cas ~ log10(p.Cymene+1) + (1|gr(Species, cov = A)))
per.p <- bf(per ~ log10(p.Cymene+1) + (1|gr(Species, cov = A)))
sow.p <- bf(sow ~ log10(p.Cymene+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.cym.p <- brm(cym.p +
        	cas.p +
        	per.p +
        	sow.p +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##beta.Caryophyllene
phy.c <- bf(log10(Caryophyllene+1) ~ (1|gr(Species, cov = A)))

##species representation
cas.c <- bf(cas ~ log10(Caryophyllene+1) + (1|gr(Species, cov = A)))
per.c <- bf(per ~ log10(Caryophyllene+1) + (1|gr(Species, cov = A)))
sow.c <- bf(sow ~ log10(Caryophyllene+1) + (1|gr(Species, cov = A)))

##run all models together at a time
m.car.b <- brm(phy.c +
        	cas.c +
        	per.c +
        	sow.c +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = dat2,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##model for MDS2
mds.2 <- bf(MDS2 ~ (1|gr(Species, cov = A)))

##species representation
cas.2 <- bf(cas ~ MDS2 + (1|gr(Species, cov = A)))
per.2 <- bf(per ~ MDS2 + (1|gr(Species, cov = A)))
sow.2 <- bf(sow ~ MDS2 + (1|gr(Species, cov = A)))

##run all models together at a time
m.mds.2 <- brm(mds.2 +
        	cas.2 +
        	per.2 +
        	sow.2 +
            set_rescor(FALSE), 
            family = list("gaussian",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 5000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat1,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

##deal with missing data
dat3<-dat2[c(2:12, 14:20, 22:22),c(1:1, 5:5, 260:263)]

##model for abundance
abu.s <- bf(abundance_salazar ~ (1|gr(Species, cov = As)))

##species representation
cas.s <- bf(cas ~ abundance_salazar + (1|gr(Species, cov = As)))
per.s <- bf(per ~ abundance_salazar + (1|gr(Species, cov = As)))
sow.s <- bf(sow ~ abundance_salazar + (1|gr(Species, cov = As)))

##run all models together at a time
m.abu.s <- brm(abu.s +
        	cas.s +
        	per.s +
        	sow.s +
            set_rescor(FALSE), 
            family = list("poisson",
            zero_inflated_beta(),
            zero_inflated_beta(),
            zero_inflated_beta()), 
                  data2 = list(As = As),
                  iter  = 100000, warmup = 50000, thin = 200, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = dat3,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_v3.RData")

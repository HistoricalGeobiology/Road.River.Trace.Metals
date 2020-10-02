library(dplyr)
library(rpart)
library(rpart.plot)

# Load data file
load("Sperling.Road.River.Data.RData")

# Calculations for authigenic enrichments (Rudnick & Gao for Mo and U)
crustal.Al <- (15.40/(26.981539*2+15.999*3))*(26.981539*2) # calculation from Al2O3 wt% to Al wt%
crustal.Mo <- 1.1
crustal.U <- 2.7
geochem$Mo.det <- (crustal.Mo/crustal.Al) * geochem$Al..wt..
geochem$Mo.auth <- geochem$Mo..ppm. - geochem$Mo.det 
geochem$U.det <- (crustal.U/crustal.Al) * geochem$Al..wt..
geochem$U.auth <- geochem$U..ppm. - geochem$U.det 
geochem$Mo_U.auth <- geochem$Mo.auth/geochem$U.auth 

#Add Mo/TOC and U/TOC values
geochem$Mo.TOC <- geochem$Mo..ppm. / geochem$TOC..wt..
geochem$U.TOC <- geochem$U..ppm. / geochem$TOC..wt..

#Add FeP/FeHR ratio where all oxides are considered to be derived from pyrite
geochem$OxideAndPyrite <- ((geochem$Fe.py..wt.. + geochem$Fe.ox..wt..)/geochem$FeHR)

#Filter by anoxic samples, then separate into euxinic and ferruginous based on iron speciation and Mo/U ratios
geochem.anox <- filter(geochem, FeHR.FeT >= 0.38 | FeT.Al >= 0.53 & !(is.na(interpreted_age)))
geochem.anox$euxinicFe[geochem.anox$Fe.py.FeHR >= 0.7] <- 1
geochem.anox$euxinicFe[geochem.anox$Fe.py.FeHR < 0.7] <- 0

geochem.anox$euxinicOxideAndPyrite[geochem.anox$OxideAndPyrite >= 0.7] <- 1
geochem.anox$euxinicOxideAndPyrite[geochem.anox$OxideAndPyrite < 0.7] <- 0

geochem.anox$euxinicMoU[geochem.anox$Mo_U.auth >= 4] <- 1
geochem.anox$euxinicMoU[geochem.anox$Mo_U.auth < 4] <- 0

#Separate into time bins for statistical analyses
Pre467 <-filter(geochem.anox, interpreted_age>=467)
Middle <-filter(geochem.anox, interpreted_age<467 & interpreted_age>408)
Post408 <-filter(geochem.anox, interpreted_age<=408)

#Create data frame just containing Peel River samples and samples for RI-07-07A drillcore
Peel <- filter(geochem.anox, !(is.na(Figure.2.Composite.Height)))
PeelPre467 <- filter(Pre467, !(is.na(Figure.2.Composite.Height)))
PeelMiddle <- filter(Middle, !(is.na(Figure.2.Composite.Height)))
PeelPost408 <- filter(Post408, !(is.na(Figure.2.Composite.Height)))

###########
##Global analyses
###########

#Calculate proportions of euxinic samples in the global dataset and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
Pre467FeSpeciationChi <- c(203,82)
MiddleFeSpeciationChi <- c(586,174)
dat <- rbind(Pre467FeSpeciationChi, MiddleFeSpeciationChi)
chisq.test(dat, correct=FALSE)

MiddleFeSpeciationChi <- c(586,174)
Post408FeSpeciationChi <- c(92,130)
dat <- rbind(MiddleFeSpeciationChi, Post408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

Pre467FeSpeciationChi <- c(203,82)
Post408FeSpeciationChi <- c(92,130)
dat <- rbind(Pre467FeSpeciationChi, Post408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

#Decision tree for best split on proportion euxinic with iron speciation in the global dataset
DecTree <- rpart(geochem.anox$euxinicFe ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Euxinic-Ferruginous iron speciation classification tree")


#Calculate proportions of euxinic samples in the global dataset (for with oxides considered to represent weathered pyrite) and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
#For all analyses, determine real p-value if R shortens to <2.2E-16
Pre467FeSpeciationOxidesChi <- c(92,175)
MiddleFeSpeciationOxidesChi <- c(352,275)
dat <- rbind(Pre467FeSpeciationOxidesChi, MiddleFeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

MiddleFeSpeciationOxidesChi <- c(352,275)
Post408FeSpeciationOxidesChi <- c(30,143)
dat <- rbind(MiddleFeSpeciationOxidesChi, Post408FeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

Pre467FeSpeciationOxidesChi <- c(92,175)
Post408FeSpeciationOxidesChi <- c(30,143)
dat <- rbind(Pre467FeSpeciationOxidesChi, Post408FeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(geochem.anox$euxinicOxideAndPyrite ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Euxinic-Ferruginous iron speciation (all oxides = pyrite) classification tree")

#Calculate proportions of euxinic samples in the global dataset using Mo/U ratios and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
Pre467euxinicMoUChi <- c(272,13)
MiddleeuxinicMoUChi <- c(557,57)
dat <- rbind(Pre467euxinicMoUChi, MiddleeuxinicMoUChi)
chisq.test(dat, correct=FALSE)

MiddleeuxinicMoUChi <- c(557,57)
Post408euxinicMoUchi <- c(200,125)
dat <- rbind(MiddleeuxinicMoUChi, Post408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

Pre467euxinicMoUChi <- c(272,13)
Post408euxinicMoUchi <- c(200,125)
dat <- rbind(Pre467euxinicMoUChi, Post408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(geochem.anox$euxinicMoU ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Euxinic-Ferruginous Mo/U classification tree")

############
##Peel River analyses
############

#Calculate proportions of euxinic samples on the Peel River and RI-07-07A and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
PeelPre467FeSpeciationChi <- c(173,34)
PeelMiddleFeSpeciationChi <- c(317,47)
dat <- rbind(PeelPre467FeSpeciationChi, PeelMiddleFeSpeciationChi)
chisq.test(dat, correct=FALSE)

PeelMiddleFeSpeciationChi <- c(317,47)
PeelPost408FeSpeciationChi <- c(70,54)
dat <- rbind(PeelMiddleFeSpeciationChi, PeelPost408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

PeelPre467FeSpeciationChi <- c(173,34)
PeelPost408FeSpeciationChi <- c(70,54)
dat <- rbind(PeelPre467FeSpeciationChi, PeelPost408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

#Decision tree for best split on proportion euxinic with iron speciation on the Peel River
DecTree <- rpart(Peel$euxinicFe ~ interpreted_age, method="class", data=Peel)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Peel River Euxinic-Ferruginous iron speciation classification tree")

#Calculate proportions of euxinic samples on the Peel River and RI-07-07A (with oxides considered to represent weathered pyrite) and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
PeelPre467FeSpeciationOxidesChi <- c(80,127)
PeelMiddleFeSpeciationOxidesChi <- c(215,149)
dat <- rbind(PeelPre467FeSpeciationOxidesChi, PeelMiddleFeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

PeelMiddleFeSpeciationOxidesChi <- c(215,149)
PeelPost408FeSpeciationChi <- c(124,105)
dat <- rbind(PeelMiddleFeSpeciationChi, PeelPost408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

PeelPre467FeSpeciationOxidesChi <- c(80,127)
PeelPost408FeSpeciationChi <- c(124,105)
dat <- rbind(PeelPre467FeSpeciationChi, PeelPost408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

#Decision tree for best split on proportion euxinic samples from the Peel River with iron speciation (all oxides considered to represent weathered pyrite) on the Peel River
DecTree <- rpart(Peel$euxinicOxideAndPyrite ~ interpreted_age, method="class", data=Peel)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Peel River Euxinic-Ferruginous iron speciation classification tree (iron = pyrite)")

#Calculate proportions of euxinic samples on the Peel River using Mo/U ratios and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
PeelPre467euxinicMoUChi <- c(222,9)
PeelMiddleeuxinicMoUChi <- c(326,38)
dat <- rbind(PeelPre467euxinicMoUChi, PeelMiddleeuxinicMoUChi)
chisq.test(dat, correct=FALSE)

PeelMiddleeuxinicMoUChi <- c(326,38)
PeelPost408euxinicMoUchi <- c(59,65)
dat <- rbind(PeelMiddleeuxinicMoUChi, PeelPost408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

PeelPre467euxinicMoUChi <- c(222,9)
PeelPost408euxinicMoUchi <- c(59,65)
dat <- rbind(PeelPre467euxinicMoUChi, PeelPost408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

#Decision tree for best split on proportion euxinic with Mo/U ratios
DecTree <- rpart(Peel$euxinicMoU ~ interpreted_age, method="class", data=Peel)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Peel River Euxinic-Ferruginous Mo/U classification tree")



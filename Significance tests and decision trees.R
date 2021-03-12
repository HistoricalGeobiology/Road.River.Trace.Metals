library(dplyr)
library(rpart)
library(rpart.plot)

# Load data file
load("Sperling.Road.River.Data.2021.RData")


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

# Binary coding of whether an anoxic sample is dominated by detrital P or authigenic P, based on a cutoff of 1000 pppm / 0.1 weight percent P
geochem.anox$authigenic.P[geochem.anox$P..ppm. >= 1000] <- 1
geochem.anox$authigenic.P[geochem.anox$P..ppm. < 1000] <- 0


#Add P/Al and Corg/Ptotal to geochem.anox, eqn from Algeo and Li 2020 GCA, multiplied by 10,000 to balance units
geochem.anox$Corg_P <- ((geochem.anox$TOC..wt../12)/(geochem.anox$P..ppm./30.97))*10000
geochem.anox$P_Al <- (geochem.anox$P..ppm./geochem.anox$Al..wt..)

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
#(no, yes)
Pre467FeSpeciationChi <- c(202,82)
MiddleFeSpeciationChi <- c(776,256)
dat <- rbind(Pre467FeSpeciationChi, MiddleFeSpeciationChi)
chisq.test(dat, correct=FALSE)

MiddleFeSpeciationChi <- c(776,256)
Post408FeSpeciationChi <- c(106,198)
dat <- rbind(MiddleFeSpeciationChi, Post408FeSpeciationChi)
chisq.test(dat, correct=FALSE)

Pre467FeSpeciationChi <- c(202,82)
Post408FeSpeciationChi <- c(106,198)
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
MiddleFeSpeciationOxidesChi <- c(422,366)
dat <- rbind(Pre467FeSpeciationOxidesChi, MiddleFeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

MiddleFeSpeciationOxidesChi <- c(422,366)
Post408FeSpeciationOxidesChi <- c(41,214)
dat <- rbind(MiddleFeSpeciationOxidesChi, Post408FeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

Pre467FeSpeciationOxidesChi <- c(92,175)
Post408FeSpeciationOxidesChi <- c(41,214)
dat <- rbind(Pre467FeSpeciationOxidesChi, Post408FeSpeciationOxidesChi)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(geochem.anox$euxinicOxideAndPyrite ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Euxinic-Ferruginous iron speciation (all oxides = pyrite) classification tree")

#Calculate proportions of euxinic samples in the global dataset using Mo/U ratios and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
Pre467euxinicMoUChi <- c(302,10)
MiddleeuxinicMoUChi <- c(778,102)
dat <- rbind(Pre467euxinicMoUChi, MiddleeuxinicMoUChi)
chisq.test(dat, correct=FALSE)

MiddleeuxinicMoUChi <- c(778,102)
Post408euxinicMoUchi <- c(520,279)
dat <- rbind(MiddleeuxinicMoUChi, Post408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

Pre467euxinicMoUChi <- c(302,10)
Post408euxinicMoUchi <- c(520,279)
dat <- rbind(Pre467euxinicMoUChi, Post408euxinicMoUchi)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(geochem.anox$euxinicMoU ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Euxinic-Ferruginous Mo/U classification tree")

#Test differences in Corg/Ptotal

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(Pre467$Corg_P)
shapiro.test(Middle$Corg_P)
shapiro.test(Post408$Corg_P)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(Pre467$Corg_P, Middle$Corg_P)
wilcox.test(Middle$Corg_P,Post408$Corg_P)
wilcox.test(Pre467$Corg_P, Post408$Corg_P)

#Regression tree for Corg/Ptotal in global dataset
RegTreeCorgPtotal <- rpart(geochem.anox$Corg_P ~ interpreted_age, method="anova", data=geochem.anox)
printcp(RegTreeCorgPtotal)
summary(RegTreeCorgPtotal)
rpart.plot(RegTreeCorgPtotal, digits=-4, main="Global Corg/Ptotal (anoxic samples) regression tree")


#Calculate proportions of anoxic samples where P is dominated by detrital P (<1000 ppm P) versus samples where P is dominated by authigenic P (>1000 ppm P), and conduct Chi square tests
#(no,yes)
Pre467anoxicPauth <- c(178,147)
MiddleanoxicPauth <- c(829,131)
dat <- rbind(Pre467anoxicPauth, MiddleanoxicPauth)
chisq.test(dat, correct=FALSE)

MiddleanoxicPauth <- c(829,131)
Post408anoxicPauth <- c(996,55)
dat <- rbind(MiddleanoxicPauth, Post408anoxicPauth)
chisq.test(dat, correct=FALSE)

Pre467anoxicPauth <- c(178,147)
Post408anoxicPauth <- c(996,55)
dat <- rbind(Pre467anoxicPauth, Post408anoxicPauth)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(geochem.anox$authigenic.P ~ interpreted_age, method="class", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Detrital-Authigenic P classification tree")

###To cut
DecTree <- rpart(geochem.anox$P..ppm. ~ interpreted_age, method="anova", data=geochem.anox)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Detrital-Authigenic P classification tree")


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

#Test differences in Corg/Ptotal

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PeelPre467$Corg_P)
shapiro.test(PeelMiddle$Corg_P)
shapiro.test(PeelPost408$Corg_P)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(PeelPre467$Corg_P, PeelMiddle$Corg_P)
wilcox.test(PeelMiddle$Corg_P,PeelPost408$Corg_P)
wilcox.test(PeelPre467$Corg_P, PeelPost408$Corg_P)

#Regression tree for Corg/Ptotal in Peel River section
RegTreeCorgPtotal <- rpart(Peel$Corg_P ~ interpreted_age, method="anova", data=Peel)
printcp(RegTreeCorgPtotal)
summary(RegTreeCorgPtotal)
rpart.plot(RegTreeCorgPtotal, digits=-4, main="Peel River Corg/Ptotal (anoxic samples) regression tree")


#Calculate proportions of anoxic samples where P is dominated by detrital P (<1000 ppm P) versus samples where P is dominated by authigenic P (>1000 ppm P), and conduct Chi square tests
#(no,yes)
PeelPre467anoxicPauth <- c(117,114)
PeelMiddleanoxicPauth <- c(277,84)
dat <- rbind(PeelPre467anoxicPauth, PeelMiddleanoxicPauth)
chisq.test(dat, correct=FALSE)

PeelMiddleanoxicPauth <- c(277,84)
PeelPost408anoxicPauth <- c(118,6)
dat <- rbind(PeelMiddleanoxicPauth, PeelPost408anoxicPauth)
chisq.test(dat, correct=FALSE)

PeelPre467anoxicPauth <- c(117,114)
PeelPost408anoxicPauth <- c(118,6)
dat <- rbind(PeelPre467anoxicPauth, PeelPost408anoxicPauth)
chisq.test(dat, correct=FALSE)

DecTree <- rpart(Peel$authigenic.P ~ interpreted_age, method="class", data=Peel)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Detrital-Authigenic P classification tree for Peel River")


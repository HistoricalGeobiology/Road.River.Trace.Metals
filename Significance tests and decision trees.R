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

#Mo analyses
#Filter for euxinic samples and samples with positive authigenic Mo
Euxinic <- filter(geochem.anox, Fe.py.FeHR >= 0.7 & Mo.auth >0)
Pre467Euxinic <- filter(Pre467, Fe.py.FeHR >= 0.7 & Mo.auth >0)
MiddleEuxinic <- filter(Middle, Fe.py.FeHR >= 0.7 & Mo.auth >0)
Post408Euxinic <- filter(Post408, Fe.py.FeHR >= 0.7 & Mo.auth >0)

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(Pre467Euxinic$Mo.auth)
shapiro.test(MiddleEuxinic$Mo.auth)
shapiro.test(Post408Euxinic$Mo.auth)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(Pre467Euxinic$Mo.auth, MiddleEuxinic$Mo.auth)
wilcox.test(MiddleEuxinic$Mo.auth,Post408Euxinic$Mo.auth)
wilcox.test(Pre467Euxinic$Mo.auth, Post408Euxinic$Mo.auth)

#Regression tree for authigenic Mo in global dataset
RegTreeMoGlobal <- rpart(Euxinic$Mo.auth ~ interpreted_age, method="anova", data=Euxinic)
printcp(RegTreeMoGlobal)
summary(RegTreeMoGlobal)
rpart.plot(RegTreeMoGlobal, digits=-4, main="Global authigenic Mo (euxinic samples) regression tree")

#Mo/TOC analyses
#No infinite Mo/TOC numbers, do not need to remove
#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(Pre467Euxinic$Mo.TOC)
shapiro.test(MiddleEuxinic$Mo.TOC)
shapiro.test(Post408Euxinic$Mo.TOC)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(Pre467Euxinic$Mo.TOC, MiddleEuxinic$Mo.TOC)
wilcox.test(MiddleEuxinic$Mo.TOC,Post408Euxinic$Mo.TOC)
wilcox.test(Pre467Euxinic$Mo.TOC, Post408Euxinic$Mo.TOC)

#Regression tree for Mo/TOC in global dataset
RegTreeMoTOCGlobal <- rpart(Euxinic$Mo.TOC ~ interpreted_age, method="anova", data=Euxinic)
printcp(RegTreeMoTOCGlobal)
summary(RegTreeMoTOCGlobal)
rpart.plot(RegTreeMoTOCGlobal, digits=-4, main="Global Mo/TOC (euxinic samples) regression tree")

#U analyses
#Filter for euxinic samples and samples with positive authigenic Mo
AnoxicU <- filter(geochem.anox, U.auth >0)
Pre467AnoxicU <- filter(Pre467, U.auth >0)
MiddleAnoxicU <- filter(Middle, U.auth >0)
Post408AnoxicU <- filter(Post408, U.auth >0)

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(Pre467AnoxicU$U.auth)
shapiro.test(MiddleAnoxicU$U.auth)
shapiro.test(Post408AnoxicU$U.auth)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(Pre467AnoxicU$U.auth, MiddleAnoxicU$U.auth)
wilcox.test(MiddleAnoxicU$U.auth,Post408AnoxicU$U.auth)
wilcox.test(Pre467AnoxicU$U.auth, Post408AnoxicU$U.auth)

#Regression tree for authigenic U in global dataset
RegTreeUGlobal <- rpart(AnoxicU$U.auth ~ interpreted_age, method="anova", data=AnoxicU)
printcp(RegTreeUGlobal)
summary(RegTreeUGlobal)
rpart.plot(RegTreeUGlobal, digits=-4, main="Global authigenic U (anoxic samples) regression tree")

#U/TOC analyses
#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(Pre467AnoxicU$U.TOC)
shapiro.test(MiddleAnoxicU$U.TOC)
shapiro.test(Post408AnoxicU$U.TOC)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(Pre467AnoxicU$U.TOC, MiddleAnoxicU$U.TOC)
wilcox.test(MiddleAnoxicU$U.TOC,Post408AnoxicU$U.TOC)
wilcox.test(Pre467AnoxicU$U.TOC, Post408AnoxicU$U.TOC)

#Regression tree for U/TOC in global dataset
RegTreeUTOCGlobal <- rpart(AnoxicU$U.TOC ~ interpreted_age, method="anova", data=AnoxicU)
printcp(RegTreeUTOCGlobal)
summary(RegTreeUTOCGlobal)
rpart.plot(RegTreeUTOCGlobal, digits=-4, main="Global U/TOC (anoxic samples) regression tree")

############
##Peel River analyses
############

#Calculate proportions of euxinic samples on the Peel River and RI-07-07A (with oxides considered to represent weathered pyrite) and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
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

#Calculate proportions of euxinic samples on the Peel River and RI-07-07A and number of euxinic versus ferruginous samples in different time bins and conduct Chi square tests
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

#Decision tree for best split on proportion euxinic with iron speciation (all oxides considered to represent weathered pyrite) on the Peel River
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

DecTree <- rpart(Peel$euxinicMoU ~ interpreted_age, method="class", data=Peel)
printcp(DecTree)
summary(DecTree)
rpart.plot(DecTree, digits=-4, main="Peel River Euxinic-Ferruginous Mo/U classification tree")

#Peel River trace metal analyses: Mo
PeelEuxinic <- filter(Peel, Fe.py.FeHR >= 0.7 & Mo.auth >0)
PeelPre467Euxinic <- filter(PeelPre467, Fe.py.FeHR >= 0.7 & Mo.auth >0)
PeelMiddleEuxinic <- filter(PeelMiddle, Fe.py.FeHR >= 0.7 & Mo.auth >0)
PeelPost408Euxinic <- filter(PeelPost408, Fe.py.FeHR >= 0.7 & Mo.auth >0)

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PeelPre467Euxinic$Mo.auth)
shapiro.test(PeelMiddleEuxinic$Mo.auth)
shapiro.test(PeelPost408Euxinic$Mo.auth)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(PeelPre467Euxinic$Mo.auth, PeelMiddleEuxinic$Mo.auth)
wilcox.test(PeelMiddleEuxinic$Mo.auth,PeelPost408Euxinic$Mo.auth)
wilcox.test(PeelPre467Euxinic$Mo.auth, PeelPost408Euxinic$Mo.auth)

#Regression tree for authigenic Mo on the Peel River
RegTreeMo <- rpart(PeelEuxinic$Mo.auth ~ interpreted_age, method="anova", data=PeelEuxinic)
printcp(RegTreeMo)
summary(RegTreeMo)
rpart.plot(RegTreeMo, digits=-4, main="Peel River authigenic Mo (euxinic samples) regression tree")

#Do not need to filter for infinite Mo/TOC values values as there are none on the Peel River
#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PeelPre467Euxinic$Mo.TOC)
shapiro.test(PeelMiddleEuxinic$Mo.TOC)
shapiro.test(PeelPost408Euxinic$Mo.TOC)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(PeelPre467Euxinic$Mo.TOC, PeelMiddleEuxinic$Mo.TOC)
wilcox.test(PeelMiddleEuxinic$Mo.TOC,PeelPost408Euxinic$Mo.TOC)
wilcox.test(PeelPre467Euxinic$Mo.TOC, PeelPost408Euxinic$Mo.TOC)

#Regression tree for Mo/TOC ratios on the Peel River
RegTreeMo.TOC <- rpart(PeelEuxinic$Mo.TOC ~ interpreted_age, method="anova", data=PeelEuxinic)
printcp(RegTreeMo.TOC)
summary(RegTreeMo.TOC)
rpart.plot(RegTreeMo.TOC, digits=-4, main="Peel River Mo/TOC (euxinic samples) regression tree")

#Peel River trace metal analyses: U
#First remove samples with U.auth <0
Peel <-filter(Peel, U.auth>0)
PeelPre467 <-filter(PeelPre467, U.auth>0)
PeelMiddle <-filter(PeelMiddle, U.auth>0)
PeelPost408 <-filter(PeelPost408, U.auth>0)

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PeelPre467$U.auth)
shapiro.test(PeelMiddle$U.auth)
shapiro.test(PeelPost408$U.auth)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(PeelPre467$U.auth, PeelMiddle$U.auth)
wilcox.test(PeelMiddle$U.auth,PeelPost408$U.auth)
wilcox.test(PeelPre467$U.auth, PeelPost408$U.auth)

#Regression tree for authigenic U on the Peel River
RegTreeU <- rpart(Peel$U.auth ~ interpreted_age, method="anova", data=Peel)
printcp(RegTreeU)
summary(RegTreeU)
rpart.plot(RegTreeU, digits=-4, main="Peel River authigenic U (anoxic samples) regression tree")

#Do not need to filter for infinite U/TOC values values as there are none on the Peel River
#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PeelPre467$U.TOC)
shapiro.test(PeelMiddle$U.TOC)
shapiro.test(PeelPost408$U.TOC)

#Normality is rejected, verified by histograms, run Mann-Whitney/Wilcoxon test
wilcox.test(PeelPre467$U.TOC, PeelMiddle$U.TOC)
wilcox.test(PeelMiddle$U.TOC,PeelPost408$U.TOC)
wilcox.test(PeelPre467$U.TOC, PeelPost408$U.TOC)

#Regression tree for Mo/TOC ratios on the Peel River
RegTreeU.TOC <- rpart(Peel$U.TOC ~ interpreted_age, method="anova", data=Peel)
printcp(RegTreeU.TOC)
summary(RegTreeU.TOC)
rpart.plot(RegTreeU.TOC, digits=-4, main="Peel River U/TOC (anoxic samples) regression tree")













#Trace metal analyses
#Filter to appropriate redox facies and time bins
#Rename original data frame
geochem_metal <-geochem

#Add Mo/TOC and U/TOC to data frame and remove infinite values
geochem_metal$Mo.TOC <- geochem_metal$Mo..ppm. / geochem_metal$TOC..wt..
geochem_metal$U.TOC <- geochem_metal$U..ppm. / geochem_metal$TOC..wt..

#Filter to appropriate redox bins (anoxic and euxinic, respectively)
geochem_metal <- filter(geochem_metal, !(is.na(interpreted_age)))
geochem.total.anox_metal <- filter(geochem_metal, FeHR.FeT > 0.38 | FeT.Al > 0.53)
euxinic_metal <-filter(geochem_metal, FeHR.FeT >= 0.38 & Fe.py.FeHR >= 0.7)

#Filter to pre- and post-408 Ma time bins
PrePragianAnoxic_metal <-filter(geochem.total.anox_metal, interpreted_age>408)
PostPragianAnoxic_metal <-filter(geochem.total.anox_metal, interpreted_age<=408)
PrePragianEuxinic_metal <-filter(euxinic_metal, interpreted_age>408)
PostPragianEuxinic_metal <-filter(euxinic_metal, interpreted_age<=408)

#Filter infinite U/TOC and Mo/TOC values
PrePragianAnoxic_metal_U.TOC <-filter(PrePragianAnoxic_metal, !(U.TOC == Inf))
PostPragianAnoxic_metal_U.TOC <-filter(PostPragianAnoxic_metal, !(U.TOC == Inf))
PrePragianEuxinic_metal_Mo.TOC <-filter(PrePragianEuxinic_metal, !(Mo.TOC == Inf))
PostPragianEuxinic_metal_Mo.TOC <-filter(PostPragianEuxinic_metal, !(Mo.TOC == Inf))

#Conduct Shapiro-Wilk tests for normality. Low values reject the null hypothesis that data come from a normal distribution
shapiro.test(PrePragianAnoxic_metal$U..ppm.)
shapiro.test(PostPragianAnoxic_metal$U..ppm.)
shapiro.test(PrePragianEuxinic_metal$Mo..ppm.)
shapiro.test(PostPragianEuxinic_metal$Mo..ppm.)
shapiro.test(PrePragianAnoxic_metal_U.TOC$U.TOC)
shapiro.test(PostPragianAnoxic_metal_U.TOC$U.TOC)
shapiro.test(PrePragianEuxinic_metal_Mo.TOC$Mo.TOC)
shapiro.test(PostPragianEuxinic_metal_Mo.TOC$Mo.TOC)

#Normality is rejected, verified by histograms run Mann-Whitney/Wilcoxon test
wilcox.test(PrePragianAnoxic_metal$U..ppm., PostPragianAnoxic_metal$U..ppm.)
wilcox.test(PrePragianEuxinic_metal$Mo..ppm., PostPragianEuxinic_metal$Mo..ppm.)

wilcox.test(PrePragianEuxinic_metal_Mo.TOC$Mo.TOC, PostPragianEuxinic_metal_Mo.TOC$Mo.TOC)
wilcox.test(PrePragianAnoxic_metal_U.TOC$U.TOC, PostPragianAnoxic_metal_U.TOC$U.TOC)

#Regression Trees for best splits in Mo, Mo/TOC, U, and U/TOC

RegTreeMo <- rpart(euxinic_metal$Mo..ppm. ~ interpreted_age, method="anova", data=euxinic_metal)
printcp(RegTreeMo)
summary(RegTreeMo)
rpart.plot(RegTreeMo, digits=-4, main="Mo regression tree")

euxinic_metal <- filter(euxinic_metal, !(Mo.TOC == Inf))
RegTreeMoTOC <- rpart(euxinic_metal$Mo.TOC ~ interpreted_age, method="anova", data=euxinic_metal)
printcp(RegTreeMoTOC)
summary(RegTreeMoTOC)
rpart.plot(RegTreeMoTOC, digits=-4, main="Mo TOC regression tree")

RegTreeU <- rpart(geochem.total.anox_metal$U..ppm. ~ interpreted_age, method="anova", data=geochem.total.anox_metal)
printcp(RegTreeU)
summary(RegTreeU)
rpart.plot(RegTreeU, digits=-4, main="U regression tree")

geochem.total.anox_metal.for.U.TOC <- filter(geochem.total.anox_metal, !(U.TOC == Inf))
RegTreeUTOC <- rpart(U.TOC ~ interpreted_age, method="anova", data=geochem.total.anox_metal.for.U.TOC)
printcp(RegTreeUTOC)
summary(RegTreeUTOC)
rpart.plot(RegTreeUTOC, digits=-4, main="U TOC regression tree")


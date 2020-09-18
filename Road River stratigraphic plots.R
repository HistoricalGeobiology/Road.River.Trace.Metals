library(readr)
library(dplyr)
load("Sperling.Road.River.Data.RData")

#Main plot (Fig. 2 of main text)
data <- filter(geochem, !(is.na(Figure.2.Composite.Height)))

#Denote whether samples have FeHR/FeT ratios >1
data$Fe.above1[data$FeHR.FeT > 1] <- 1
data$Fe.above1[data$FeHR.FeT <= 1]  <- 0

# Calculations for authigenic enrichments (Rudnick & Gao 2003 for Mo and U)
crustal.Al <- (15.40/(26.981539*2+15.999*3))*(26.981539*2) # calculation from Al2O3 wt% to Al wt%
crustal.Mo <- 1.1
crustal.U <- 2.7
crustal.Cr <- 92
crustal.Ni <- 47
crustal.V <- 97
crustal.Zn <- 67

#Define authigenic and detrital values
data$Mo.det <- (crustal.Mo/crustal.Al) * data$Al..wt..
data$Mo.auth <- data$Mo..ppm. - data$Mo.det 
data$U.det <- (crustal.U/crustal.Al) * data$Al..wt..
data$U.auth <- data$U..ppm. - data$U.det 
data$Mo_U.auth <- data$Mo.auth/data$U.auth 
data$Cr.det <- (crustal.Cr/crustal.Al) * data$Al..wt..
data$Cr.auth <- data$Cr..ppm. - data$Cr.det 
data$Ni.det <- (crustal.Ni/crustal.Al) * data$Al..wt..
data$Ni.auth <- data$Ni..ppm. - data$Ni.det 
data$V.det <- (crustal.V/crustal.Al) * data$Al..wt..
data$V.auth <- data$V..ppm. - data$V.det 
data$Zn.det <- (crustal.Zn/crustal.Al) * data$Al..wt..
data$Zn.auth <- data$Zn..ppm. - data$Zn.det 

#Add FeP/FeHR ratio where all oxides are considered to be derived from pyrite
data$OxideAndPyrite <- ((data$Fe.py..wt.. + data$Fe.ox..wt..)/data$FeHR)


#Plot Figure 2 in Main text. In iron speciation plots, points with FeHR/FeT > 1 or FeT < 0.5 weight percent are plotted as open circles
quartz(w=15, h=12)
layout(matrix(c(1:6),1,6))
plot(data$FeHR.FeT, data$Figure.2.Composite.Height, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylab="", yaxt='n', ylim=c(1096.2, 2607.84), xlim=c(0,1.0), abline(v=0.38), cex.lab=1.5, xlab="FeHR/FeT")
plot(data$Fe.py.FeHR, data$Figure.2.Composite.Height, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(1096.2, 2607.84), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(data$TOC, data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(data$Mo.auth, data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(data$U.auth, data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(data$Mo_U.auth, data$Figure.2.Composite.Height, yaxs="i", pch=19,  ylim=c(1096.2, 2607.84), xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Then scale 80%, 80%, 80%, 97% and 99% in illustrator to match stratigraphic column

#Supplemental plots

#Add metal/TOC and EFs (using Turekian and Wedepohl, 1961 average shale values) to data frame and remove infinite values
data$Mo.TOC <- data$Mo..ppm. / data$TOC..wt..
data$U.TOC <- data$U..ppm. / data$TOC..wt..


#Supplemental plot 2: Major element data, metal/TOC data, and FeP/FeHR ratios with all iron oxides considered to represent oxidized pyrite. 
quartz(w=15, h=12)
layout(matrix(c(1:6),1,6))
plot(data$OxideAndPyrite, data$Figure.2.Composite.Height, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(1096.2, 2607.84), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(data$Si..wt.., data$Figure.2.Composite.Height, yaxs="i", pch=19, ylab="",yaxt='n', ylim=c(1096.2, 2607.84), xlim=c(0,50), cex.lab=1.5, xlab="Si wt %")
plot(data$Ca..wt..,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,35), cex.lab=1.5, xlab="Ca wt %", ylab="",yaxt='n')
plot(data$Al..wt..,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,12),cex.lab=1.5, xlab="Al wt %", ylab="",yaxt='n')
plot(data$Mo.TOC,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,75), cex.lab=1.5, xlab="Mo/TOC", ylab="",yaxt='n')
plot(data$U.TOC,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,50), cex.lab=1.5, xlab="U/TOC", ylab="",yaxt='n')

#Supplemental plot 3: Additional redox-sensitive elements
quartz(w=15, h=12)
layout(matrix(c(1:6),1,6))
plot(data$Cr..ppm., data$Figure.2.Composite.Height, yaxs="i", pch=19, ylab="",yaxt='n', ylim=c(1096.2, 2607.84), xlim=c(0, 250), cex.lab=1.5, xlab="Cr (ppm)")
plot(data$Ni..ppm., data$Figure.2.Composite.Height, yaxs="i", pch=19, ylab="",yaxt='n', ylim=c(1096.2, 2607.84), xlim=c(0,250), cex.lab=1.5, xlab="Ni (ppm)")
plot(data$V..ppm.,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,4000), cex.lab=1.5, xlab="V (ppm)", ylab="",yaxt='n')
plot(data$Zn..ppm.,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,5000),cex.lab=1.5, xlab="Zn (ppm)", ylab="",yaxt='n')
plot(data$P..ppm.,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,20000), cex.lab=1.5, xlab="P (ppm)", ylab="",yaxt='n')
plot(data$Mn..ppm.,data$Figure.2.Composite.Height, yaxs="i", pch=19, ylim=c(1096.2, 2607.84), xlim=c(0,1000), cex.lab=1.5, xlab="Mn (ppm)", ylab="",yaxt='n')

#Supplemental Plot: J1728
dataJ1728 <-filter(data, section_name=="J1728")
quartz(w=15, h=9.85)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1728$delta_13C.org..permil., dataJ1728$height.depth.in.meters, pch=19, ylim=c(0, 306.2), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1728$FeHR.FeT, dataJ1728$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 306.2), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1728$Fe.py.FeHR, dataJ1728$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 306.2), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1728$TOC, dataJ1728$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 306.2), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1728$Mo.auth, dataJ1728$height.depth.in.meters, yaxs="i", ylim=c(0, 306.2), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1728$U.auth, dataJ1728$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 306.2),  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1728$Mo_U.auth, dataJ1728$height.depth.in.meters, yaxs="i", ylim=c(0, 306.2), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')


#Supplemental Plot: TF-17-03
data17TF03 <-filter(data, section_name=="17-TF-03")
quartz(w=15, h=5.53)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(data17TF03$delta_13C.org..permil., data17TF03$height.depth.in.meters, pch=19, ylim=c(0, 162), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(data17TF03$FeHR.FeT, data17TF03$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 162), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(data17TF03$Fe.py.FeHR, data17TF03$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 162), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(data17TF03$TOC, data17TF03$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 162), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(data17TF03$Mo.auth, data17TF03$height.depth.in.meters, yaxs="i", ylim=c(0, 162), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(data17TF03$U.auth, data17TF03$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 162),  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(data17TF03$Mo_U.auth, data17TF03$height.depth.in.meters, yaxs="i", ylim=c(0, 162), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')


#Supplemental Plot: J1727
dataJ1727 <-filter(data, section_name=="J1727")
quartz(w=15, h=8.079)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1727$delta_13C.org..permil., dataJ1727$height.depth.in.meters, pch=19, ylim=c(0, 247.6), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1727$FeHR.FeT, dataJ1727$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 247.6), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1727$Fe.py.FeHR, dataJ1727$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 247.6), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1727$TOC, dataJ1727$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 247.6), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1727$Mo.auth, dataJ1727$height.depth.in.meters, yaxs="i", ylim=c(0, 247.6), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1727$U.auth, dataJ1727$height.depth.in.meters, yaxs="i", ylim=c(0, 247.6), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1727$Mo_U.auth, dataJ1727$height.depth.in.meters, yaxs="i", ylim=c(0, 247.6), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Supplemental Plot: J1611
dataJ1611 <-filter(data, section_name=="J1611")
quartz(w=15, h=1.396)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1611$delta_13C.org..permil., dataJ1611$height.depth.in.meters, pch=19, ylim=c(0, 23.9), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1611$FeHR.FeT, dataJ1611$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 23.9), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1611$Fe.py.FeHR, dataJ1611$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 23.9), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1611$TOC, dataJ1611$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 23.9), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1611$Mo.auth, dataJ1611$height.depth.in.meters, yaxs="i", ylim=c(0, 23.9), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1611$U.auth, dataJ1611$height.depth.in.meters, yaxs="i", ylim=c(0, 23.9), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1611$Mo_U.auth, dataJ1611$height.depth.in.meters, yaxs="i", ylim=c(0, 23.9), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Supplemental Plot: J1729
dataJ1729 <-filter(data, section_name=="J1729")
quartz(w=15, h=1.487)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1729$delta_13C.org..permil., dataJ1729$height.depth.in.meters, pch=19, ylim=c(0, 26.9), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1729$FeHR.FeT, dataJ1729$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 26.9), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1729$Fe.py.FeHR, dataJ1729$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 26.9), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1729$TOC, dataJ1729$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 26.9), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1729$Mo.auth, dataJ1729$height.depth.in.meters, yaxs="i", ylim=c(0, 26.9), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1729$U.auth, dataJ1729$height.depth.in.meters, yaxs="i", ylim=c(0, 26.9), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1729$Mo_U.auth, dataJ1729$height.depth.in.meters, yaxs="i", ylim=c(0, 26.9), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Supplemental Plot: J1518
dataJ1518 <-filter(data, section_name=="J1518")
quartz(w=15, h=9.32)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1518$delta_13C.org..permil., dataJ1518$height.depth.in.meters, pch=19, ylim=c(0, 288.7), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1518$FeHR.FeT, dataJ1518$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 288.7), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1518$Fe.py.FeHR, dataJ1518$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 288.7), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1518$TOC, dataJ1518$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 288.7), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1518$Mo.auth, dataJ1518$height.depth.in.meters, yaxs="i", ylim=c(0, 288.7), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1518$U.auth, dataJ1518$height.depth.in.meters, yaxs="i", ylim=c(0, 288.7), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1518$Mo_U.auth, dataJ1518$height.depth.in.meters, yaxs="i", ylim=c(0, 288.7), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Supplemental Plot: 15-TF-05 and 15-TF-07
data15TF0507 <-filter(data, section_name=="15-TF-05" | section_name=="15-TF-07")
#Make composite stratigraphic heights for section
data15TF0507$comp.height= ifelse(data15TF0507$section_name=="15-TF-07", data15TF0507$height.depth.in.meters, data15TF0507$height.depth.in.meters+44.5)  

quartz(w=15, h=7.865)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(data15TF0507$delta_13C.org..permil., data15TF0507$comp.height, pch=19, ylim=c(-.1, 239.5), yaxs="i", ylab="", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(data15TF0507$FeHR.FeT, data15TF0507$comp.height, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(-.1, 239.5), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(data15TF0507$Fe.py.FeHR, data15TF0507$comp.height, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(-.1, 239.5), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(data15TF0507$TOC, data15TF0507$comp.height, yaxs="i", pch=19, ylim=c(0, 239.5), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(data15TF0507$Mo.auth, data15TF0507$comp.height, yaxs="i", ylim=c(-.1, 239.5), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(data15TF0507$U.auth, data15TF0507$comp.height, yaxs="i", ylim=c(-.1, 239.5), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(data15TF0507$Mo_U.auth, data15TF0507$comp.height, yaxs="i", ylim=c(-.1, 239.5), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')

#Supplemental Plot: J1609
dataJ1609 <-filter(data, section_name=="J1609")
quartz(w=15, h=3.482)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1609$delta_13C.org..permil., dataJ1609$height.depth.in.meters, pch=19, ylim=c(0, 92.9), ylab="", yaxs="i", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1609$FeHR.FeT, dataJ1609$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 92.9), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1609$Fe.py.FeHR, dataJ1609$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 92.9), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1609$TOC, dataJ1609$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 92.9), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1609$Mo.auth, dataJ1609$height.depth.in.meters, yaxs="i", ylim=c(0, 92.9), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1609$U.auth, dataJ1609$height.depth.in.meters, yaxs="i", ylim=c(0, 92.9), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1609$Mo_U.auth, dataJ1609$height.depth.in.meters, yaxs="i", ylim=c(0, 92.9), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')


#Supplemental Plot: J1610
dataJ1610 <-filter(data, section_name=="J1610")
quartz(w=15, h=9.545)
#Height should be height of strat column in illustrator + 0.67 inches
layout(matrix(c(1:7),1,7))
par(mai=c(.67,.2,0,.2))
plot(dataJ1610$delta_13C.org..permil., dataJ1610$height.depth.in.meters, pch=19, ylim=c(0, 295.6), yaxs="i", ylab="", yaxt='n',xlim=c(-32.5, -25.5), cex.lab=1.5, xlab="Delta 13C organic")
plot(dataJ1610$FeHR.FeT, dataJ1610$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), ylim=c(0, 295.6), ylab="", yaxt='n', abline(v=0.38), xlim=c(0,1.0), cex.lab=1.5,  xlab="FeHR/FeT")
plot(dataJ1610$Fe.py.FeHR, dataJ1610$height.depth.in.meters, yaxs="i", pch= ifelse((data$FeT<=0.501 | data$Fe.above1==1), 1, 19), xlim=c(0,1.0), abline(v=0.7), cex.lab=1.5, xlab="FeP/FeHR", ylab="",yaxt='n')
plot(dataJ1610$TOC, dataJ1610$height.depth.in.meters, yaxs="i", pch=19, ylim=c(0, 295.6), cex.lab=1.5, xlab="TOC (wt%)", ylab="",yaxt='n')
plot(dataJ1610$Mo.auth, dataJ1610$height.depth.in.meters, yaxs="i", ylim=c(0, 295.6), pch=19,  cex.lab=1.5, xlab="Mo auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1610$U.auth, dataJ1610$height.depth.in.meters, yaxs="i", ylim=c(0, 295.6), pch=19,  cex.lab=1.5, xlab="U auth. (ppm)", ylab="",yaxt='n')
plot(dataJ1610$Mo_U.auth, dataJ1610$height.depth.in.meters, yaxs="i", ylim=c(0, 295.6), pch=19, xlim=c(0,20), abline(v=4), cex.lab=1.5, xlab="Mo auth. / U auth.", ylab="",yaxt='n')





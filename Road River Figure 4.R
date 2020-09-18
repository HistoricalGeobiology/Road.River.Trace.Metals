library(ggplot2)
library(deeptime)
library(dplyr)
library(binom)
library(simpleboot)

# Set working directory to wherever you have stored the invweight.R function! Plots will also be saved here. 
setwd()

# Set age parameters for whole analysis 
bin_size <- 10
start_age <- 490
end_age <- 360

# Rename abbreviated stage names in deeptime to fit intended box size
stages.renamed <- deeptime::stages
stages.renamed$abbr[69:91] <- c("Fm", "Frs", "Gv", "E", "Em", "P", "Lc", "P", "L", "", "H", "S", "T", "A", "R", "", "Kt", "S", "Drr", "D", "Fl", "Trm", "10")

# Import data 
load("Sperling.Road.River.Data.RData")

# Filter main dataset to generate anoxic and euxinic dataframes for plotting. 
geochem <- filter(geochem, !(is.na(interpreted_age)))
geochem.anox <- filter(geochem, FeHR.FeT >= 0.38 | FeT.Al >= 0.53)

geochem.anox$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2
geochem.anox <- filter(geochem.anox, !(is.na(time.bin)))
geochem.anox$U.TOC <- geochem.anox$U..ppm./geochem.anox$TOC..wt..

geochem.eux <- filter(geochem.anox, Fe.py.FeHR >= 0.7)
geochem.eux$Mo.TOC <- geochem.eux$Mo..ppm./geochem.eux$TOC..wt..

# Calculations for authegenic enrichments (Rudnick & Gao for Mo and U)
crustal.Al <- (15.40/(26.981539*2+15.999*3))*(26.981539*2) # calculation from Al2O3 wt% to Al wt%
crustal.Mo <- 1.1
crustal.U <- 2.7

geochem.anox$U.det <- (crustal.U/crustal.Al) * geochem.anox$Al..wt..
geochem.anox$U.auth <- geochem.anox$U..ppm. - geochem.anox$U.det 

geochem.eux$Mo.det <- (crustal.Mo/crustal.Al) * geochem.eux$Al..wt..
geochem.eux$Mo.auth <- geochem.eux$Mo..ppm. - geochem.eux$Mo.det 

# Generate separate ggplot objects (titles explanatory)

U.plot <- ggplot(geochem.anox, aes(x=time.bin, y=U.auth, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4) +
  geom_boxplot(width=5)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-2,65),expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression("U"[auth]*" (ppm)"))+xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

U.TOC.plot <- ggplot(filter(geochem.anox, !(is.na(U.TOC) | U.TOC ==Inf)), aes(x=time.bin, y=U.TOC, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4) +
  geom_boxplot(width=5)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-1,32),expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("U/TOC")+xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())


Mo.plot <- ggplot(geochem.eux, aes(x=time.bin, y=Mo.auth, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4 ) +
  geom_boxplot(width=5)+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(-3,155),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size = list(5, 8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression("Mo"[auth]*" (ppm)"))+xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

Mo.TOC.plot <- ggplot(filter(geochem.eux, !(is.na(Mo.TOC) | Mo.TOC ==Inf)), aes(x=time.bin, y=Mo.TOC, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4 ) +
  geom_boxplot(width=5)+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(-1,32),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size = list(5, 8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("Mo/TOC")+xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Combine to make version of Fig 4 with just boxplots (no invweight analysis of global data)
fig4.noboot <- ggarrange2(U.plot, U.TOC.plot, Mo.plot, Mo.TOC.plot, ncol=2, heights = c(1,1))

####################################################
## Spatial bootstrap of global datasets
####################################################

# Source invweight function - same function as in Keller & Schoene (Nature, 2012)
source("invweight.R")

####################################################
## Authigenic URANIUM (PPM)
####################################################

## Filter for anoxic U only
U.anox <-  geochem.anox  %>%
           filter(!is.na(U.auth))

## Geospatial weighting coefficients for U in anoxic shale [function "invweight" uses declustering eqn from Keller et al. 2012]
U.anox$k <- invweight(lat = U.anox$lat_dec,
                      lon = U.anox$long_dec,
                      age = U.anox$interpreted_age)

## Visualize geographic weighting
ggplot(U.anox, aes(x=long_dec, y=lat_dec))+geom_point(aes(size=k), alpha=.2)+theme_bw()

## Bin U data by age
U.anox$age.bin <-  as.numeric(cut(U.anox$interpreted_age, seq(end_age, start_age, bin_size)))

## Bootstrap anoxic Uranium
u_bootmeans <- as.numeric() 

## Loop through time bins and generate mean authigenic concentration values based on spatial bootstrap
for (bin in 1:max(U.anox$age.bin)){
  timebin <- U.anox %>%
    filter(age.bin==bin) %>%
    filter(!is.na(U.auth)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <- one.boot(timebin$U.auth, mean, R=1000, weights=1/timebin$k)
    u_bootmeans <- rbind(u_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}

## Manipulate data for plotting
uranium <- as.data.frame(cbind(u_bootmeans, (end_age-bin_size/2)+u_bootmeans[,2]*bin_size))
names(uranium) <- c("U", "bin", "time")

# Generate an individual bootstrapped mean per time bin
U.sum <- uranium %>%
  group_by(time) %>%
  summarise(U.mean = mean(U))

U.plot.w.boot <- ggplot(geochem.anox, aes(x=time.bin, y=U.auth, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4) +
  geom_boxplot(width=5)+
  annotate(geom="point", x=U.sum$time, y=U.sum$U.mean, size=5, fill = "grey70", shape=21, stroke=1)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-2,65),expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression("U"[auth]*" (ppm)"))+xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
U.plot.w.boot

####################################################
## URANIUM/TOC
####################################################

##Filter for U and TOC in aanoxic shales
U_TOC.anox <- geochem.anox %>%
  filter(U..ppm. > 0) %>%
  filter(!is.na(U..ppm.)) %>%
  filter(TOC..wt.. > 0.3) %>%
  filter(!is.na(TOC..wt..))

## Geospatial weighting coefficients for U and TOC in anoxic shale [function "invweight" uses declustering eqn from Keller et al. 2012]
U_TOC.anox$k <- invweight(lat = U_TOC.anox$lat_dec,
                          lon = U_TOC.anox$long_dec,
                          age = U_TOC.anox$interpreted_age)

## Bin U_TOC data by age
U_TOC.anox$age.bin <-  as.numeric(cut(U_TOC.anox$interpreted_age, seq(end_age, start_age, bin_size)))

## Bootstrapping anoxic Uranium/TOC
u_toc_bootmeans <- as.numeric()

## Loop through time bins and generate mean U/TOC values based on spatial bootstrap
for (bin in 1:max(U_TOC.anox$age.bin)){
  timebin <- U_TOC.anox %>%
    filter(age.bin==bin) %>%
    filter(!is.na(U..ppm.)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <-one.boot(timebin$U.TOC, mean, R=1000, weights=1/timebin$k)
    u_toc_bootmeans <- rbind(u_toc_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}


## Manipulate data for plotting
uranium_toc <- as.data.frame(cbind(u_toc_bootmeans, (end_age-bin_size/2)+u_toc_bootmeans[,2]*bin_size))

names(uranium_toc) <- c("U.TOC", "bin", "time")

# Generate an individual bootstrapped mean per time bin
U.TOC.sum <- uranium_toc %>%
  group_by(time) %>%
  summarise(U.TOC.mean = mean(U.TOC))

U.TOC.plot.w.boot <- ggplot(filter(geochem.anox, !(is.na(U.TOC) | U.TOC ==Inf)), aes(x=time.bin, y=U.TOC, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4) +
  geom_boxplot(width=5)+
  annotate(geom="point", x=U.TOC.sum$time, y=U.TOC.sum$U.TOC.mean, size=5, fill = "grey70", shape=21, stroke=1)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-1,32),expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("U/TOC")+xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())


####################################################
## Authigenic MOLYBDENUM (PPM)
####################################################

##Filter for euxinic Mo only
Mo.eux <- geochem.eux %>%
  filter(!is.na(Mo.auth))

## Geospatial weighting coefficients for Mo in euxinic shale [function "invweight" uses declustering eqn from Keller et al. 2012]
Mo.eux$k <- invweight(lat = Mo.eux$lat_dec,
                      lon = Mo.eux$long_dec,
                      age = Mo.eux$interpreted_age)

ggplot(Mo.eux, aes(x=long_dec, y=lat_dec))+geom_point(aes(size=k), alpha=.2)+theme_bw()

## Bin Mo data by age
Mo.eux$age.bin <-  as.numeric(cut(Mo.eux$interpreted_age, seq(end_age, start_age, bin_size)))

## Bootstrap euxinic molybdenum
Mo_bootmeans <- as.numeric() 

## Loop through time bins and generate mean authigenic concentration values based on spatial bootstrap
for (bin in 1:max(Mo.eux$age.bin)){
  timebin <- Mo.eux %>%
    filter(age.bin==bin) %>%
    filter(!is.na(Mo.auth)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <- one.boot(timebin$Mo.auth, mean, R=1000, weights=1/timebin$k)
    Mo_bootmeans <- rbind(Mo_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}

## Manipulate data for plotting
molybdenum <- as.data.frame(cbind(Mo_bootmeans, (end_age-bin_size/2)+Mo_bootmeans[,2]*bin_size))

names(molybdenum) <- c("Mo", "bin", "time")

# Generate an individual bootstrapped mean per time bin
Mo.sum <- molybdenum %>%
  group_by(time) %>%
  summarise(Mo.mean = mean(Mo))

Mo.plot.w.boot <- ggplot(geochem.eux, aes(x=time.bin, y=Mo.auth, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4 ) +
  geom_boxplot(width=5)+
  annotate(geom="point", x=Mo.sum$time, y=Mo.sum$Mo.mean, size=5, fill = "grey70", shape=21, stroke=1)+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(-3,155),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size = list(5, 8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression("Mo"[auth]*" (ppm)"))+xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

####################################################
## Molybdenum/TOC
####################################################

## Filter for TOC and Mo in euxinic shales
Mo_TOC.eux <- geochem.eux %>%
  filter(Mo..ppm. > 0) %>%
  filter(!is.na(Mo..ppm.)) %>%
  filter(TOC..wt.. > 0.3) %>%
  filter(!is.na(TOC..wt..))


## Geospatial weighting coefficients for Mo and TOC in euxinic shale [function "invweight" uses declustering eqn from Keller et al. 2012]
Mo_TOC.eux$k <- invweight(lat = Mo_TOC.eux$lat_dec,
                          lon = Mo_TOC.eux$long_dec,
                          age = Mo_TOC.eux$interpreted_age)

## Add k values to Mo/TOC dataset
Mo_TOC.eux$age.bin <-  as.numeric(cut(Mo_TOC.eux$interpreted_age, seq(end_age, start_age, bin_size)))

## Bootstrapping euxinic molybdenum/TOC
Mo_toc_bootmeans <- as.numeric()

## Loop through time bins and generate mean Mo/TOC values based on spatial bootstrap
for (bin in 1:max(Mo_TOC.eux$age.bin)){
  timebin <- Mo_TOC.eux %>%
    filter(age.bin==bin) %>%
    filter(!is.na(Mo..ppm.)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <-one.boot(timebin$Mo.TOC, mean, R=1000, weights=1/timebin$k)
    Mo_toc_bootmeans <- rbind(Mo_toc_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}


## Manipulate data for plotting
molybdenum_toc <- as.data.frame(cbind(Mo_toc_bootmeans, (end_age-bin_size/2)+Mo_toc_bootmeans[,2]*bin_size))

names(molybdenum_toc) <- c("Mo.TOC", "bin", "time")

# Generate an individual bootstrapped mean per time bin
Mo.TOC.sum <- molybdenum_toc %>%
  group_by(time) %>%
  summarise(Mo.TOC.mean = mean(Mo.TOC))

Mo.TOC.plot.w.boot <-  ggplot(filter(geochem.eux, !(is.na(Mo.TOC) | Mo.TOC ==Inf)), aes(x=time.bin, y=Mo.TOC, group=time.bin))+
  stat_boxplot(geom ='errorbar', size=0.6, width = 4 ) +
  geom_boxplot(width=5)+
  annotate(geom="point", x=Mo.TOC.sum$time, y=Mo.TOC.sum$Mo.TOC.mean, size=5, fill = "grey70", shape=21, stroke=1)+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(-1,32),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size = list(5, 8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("Mo/TOC")+xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Road River fig 4 with bootstrap but no histograms
fig4.w.boot <- ggarrange2(U.plot.w.boot, U.TOC.plot.w.boot, Mo.plot.w.boot, Mo.TOC.plot.w.boot, ncol=2)

## ---------------------------------------------------------------------
## Generate histograms of sampling for each of the four panels in figure 4
## ---------------------------------------------------------------------

U.anox$road.riv <- NA
U.anox$road.riv[U.anox$lithostrat != "Road River Group"] <- "Global Compilation"
U.anox$road.riv[U.anox$lithostrat == "Road River Group"] <- "Road River Group"

U.plot.w.boot_sample.hist <- ggplot(U.anox, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,200), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200), labels=c("","","100", "", "200"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1), 
        axis.ticks = element_line(size=1.5), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position='none',
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
U.plot.w.boot_sample.hist

U_TOC.anox$road.riv <- NA
U_TOC.anox$road.riv[U_TOC.anox$lithostrat != "Road River Group"] <- "Global Compilation"
U_TOC.anox$road.riv[U_TOC.anox$lithostrat == "Road River Group"] <- "Road River Group"

U.TOC.plot.w.boot_sample.hist <- ggplot(U_TOC.anox, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,200), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200), labels=c("","","100", "", "200"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1), 
        axis.ticks = element_line(size=1.5), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,1.4),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
U.TOC.plot.w.boot_sample.hist

Mo.eux$road.riv <- NA
Mo.eux$road.riv[Mo.eux$lithostrat != "Road River Group"] <- "Global Compilation"
Mo.eux$road.riv[Mo.eux$lithostrat == "Road River Group"] <- "Road River Group"

Mo.plot.w.boot_sample.hist <- ggplot(Mo.eux, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,100), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200)*.5, labels=c("","","50", "", "100"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1.5), 
        axis.ticks = element_line(size=1), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position="none",
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Mo.plot.w.boot_sample.hist

Mo_TOC.eux$road.riv <- NA
Mo_TOC.eux$road.riv[Mo_TOC.eux$lithostrat != "Road River Group"] <- "Global Compilation"
Mo_TOC.eux$road.riv[Mo_TOC.eux$lithostrat == "Road River Group"] <- "Road River Group"

Mo.TOC.plot.w.boot_sample.hist <- ggplot(Mo_TOC.eux, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,100), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200)*.5, labels=c("","","50", "", "100"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(1,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1.5), 
        axis.ticks = element_line(size=1), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position="none",
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Mo.TOC.plot.w.boot_sample.hist

# Full figure 4 with bootstrap and histograms
fig4.w.boot.w.hist <- ggarrange2(U.plot.w.boot_sample.hist, U.TOC.plot.w.boot_sample.hist, U.plot.w.boot, U.TOC.plot.w.boot, Mo.plot.w.boot_sample.hist, Mo.TOC.plot.w.boot_sample.hist, Mo.plot.w.boot, Mo.TOC.plot.w.boot, ncol=2, heights=c(0.2,1,0.2,1), widths=c(1,1))

ggsave(file="Road River Figure 4.pdf", fig4.w.boot.w.hist, height=18, width=11*2)



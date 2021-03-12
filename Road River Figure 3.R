library(ggplot2)
library(dplyr)
library(deeptime)
library(egg)
library(simpleboot)

# Set working directory for plot to be saved in
setwd()

# set seed for bootstrap analyses
set.seed(408)

# Specify start and end ages
start_age <- 490
end_age <- 360
bin_size <- 10

# Adapt abbreviated stage names in deeptime package for plotting.
stages.renamed <- deeptime::stages
stages.renamed$abbr[69:91] <- c("Fm", "Frs", "Gv", "E", "Em", "P", "Lc", "P", "L", "", "H", "S", "T", "A", "R", "", "Kt", "S", "Drr", "D", "Fl", "Trm", "10")

# Source invweight function - function from Keller & Schoene (Nature, 2012)
source("invweight.R")

# Load data file
load("Sperling.Road.River.Data.2021.RData")
## ------------------------------------------------------------
# Panel A - proportion euxinic samples based on iron speciation
## ------------------------------------------------------------

# Filter to anoxic samples with ages only
geochem.anox <- filter(geochem, FeHR.FeT >= 0.38 & !(is.na(interpreted_age)) & !(is.na(lat_dec)) & !(is.na(long_dec)) &
                       interpreted_age <= start_age, interpreted_age >= end_age)


## Geospatial weighting coefficients for anoxic shale samples [function "invweight" uses declustering eqn from Keller et al. 2012]
geochem.anox$k <- invweight(lat = geochem.anox$lat_dec,
                      lon = geochem.anox$long_dec,
                      age = geochem.anox$interpreted_age)

# Divide samples into age bins
geochem.anox$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2

# Binary coding of euxinic samples
geochem.anox$euxinic.Fe[geochem.anox$Fe.py.FeHR >= 0.7] <- 1
geochem.anox$euxinic.Fe[geochem.anox$Fe.py.FeHR < 0.7] <- 0

## Bootstrap bimodal euxinia classification using invweight
eux_bootmeans <- as.numeric() 

## Loop through time bins and generate mean authigenic concentration values based on spatial bootstrap
for (bin in (seq(end_age, start_age,  bin_size)+bin_size/2)){
  timebin <- geochem.anox %>%
    filter(time.bin==bin) %>%
    filter(!is.na(euxinic.Fe)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <- one.boot(timebin$euxinic.Fe, mean, R=1000, weights=1/timebin$k)
    eux_bootmeans <- rbind(eux_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}

## Generate dataframe of results and name columns
eux.Fe <- as.data.frame(eux_bootmeans)
names(eux.Fe) <- c("Prop", "Age")

# Generate an individual bootstrapped mean per time bin
eux.sum.Fe <- eux.Fe %>%
  group_by(Age) %>%
  summarise(Prop.mean = mean(Prop), Prop.min = mean(Prop) - 2*(sd(Prop, na.rm=T)), Prop.max = mean(Prop) + 2 * (sd(Prop, na.rm=T)))

# Generate panel A plot
eux.plot.Fe <- ggplot(eux.sum.Fe, aes(x=Age, y=Prop.mean, ymin=Prop.min, ymax=Prop.max))+
  geom_pointrange(size=1.5)+
  geom_point(shape=21, size=7, stroke=1.5, fill="grey80", alpha=.99)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-0.05,1.05), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression(atop("Proportion ", paste("euxinic (Fe"[Py]*"/Fe"[HR]*")"))))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.3,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position="none",
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="text", label="A", size=13, y=0.97, x=486)
eux.plot.Fe

# Generate histogram of samples for panel A

geochem.anox$road.riv <- NA
geochem.anox$road.riv[geochem.anox$lithostrat != "Road River Group"] <- "Global Compilation"
geochem.anox$road.riv[geochem.anox$lithostrat == "Road River Group"] <- "Road River Group"

# Summarize counts to ensure histograms span full range of sample numbers
geochem.anox %>% 
  group_by(time.bin) %>%
  tally()

eux.plot.Fe_sample.hist <- ggplot(geochem.anox, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,425), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300,350,400), labels=c("","","","","200","","","", "400"))+
  ylab("n")+
  xlab("Time (Ma)")+
  guides(fill = guide_legend(keywidth = 1.8))+
  theme(plot.margin = margin(.3,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1), 
        axis.ticks = element_line(size=1.5), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=30),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(1.0,1.25),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
eux.plot.Fe_sample.hist


## ------------------------------------------------------------
# Panel B - proportion euxinic samples based on Mo/U
## ------------------------------------------------------------

# Euxinic cutoff - choose Mo/U ratio. 4 is used as default.
eux.cutoff <- 4

# Filter samples - all anoxic samples with age, Mo, U and Al data
geochem.anox_all.MoU <- filter(geochem, (FeHR.FeT >= 0.38 | FeT.Al >= 0.53 ) & !(is.na(interpreted_age)) & 
                                 !(is.na(Mo..ppm.)) &  !(is.na(U..ppm.)) & !(is.na(Al..wt..)) &  
                                 !(is.na(lat_dec)) & !(is.na(long_dec)) &
                                 interpreted_age <= start_age, interpreted_age >= end_age)


# Calculations for authegenic enrichments (Rudnick & Gao for Mo and U)
crustal.Al <- (15.40/(26.981539*2+15.999*3))*(26.981539*2) # calculation from Al2O3 wt% to Al wt%
crustal.Mo <- 1.1
crustal.U <- 2.7
geochem.anox_all.MoU$Mo.det <- (crustal.Mo/crustal.Al) * geochem.anox_all.MoU$Al..wt..
geochem.anox_all.MoU$Mo.auth <- geochem.anox_all.MoU$Mo..ppm. - geochem.anox_all.MoU$Mo.det 
geochem.anox_all.MoU$U.det <- (crustal.U/crustal.Al) * geochem.anox_all.MoU$Al..wt..
geochem.anox_all.MoU$U.auth <- geochem.anox_all.MoU$U..ppm. - geochem.anox_all.MoU$U.det 
geochem.anox_all.MoU$Mo_U.auth <- geochem.anox_all.MoU$Mo.auth/geochem.anox_all.MoU$U.auth 

## Geospatial weighting coefficients for anoxic shale samples [function "invweight" uses declustering eqn from Keller et al. 2012]
geochem.anox_all.MoU$k <- invweight(lat = geochem.anox_all.MoU$lat_dec,
                            lon = geochem.anox_all.MoU$long_dec,
                            age = geochem.anox_all.MoU$interpreted_age)

# Divide samples into age bins
geochem.anox_all.MoU$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox_all.MoU$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2

# Binary coding of euxinic samples
geochem.anox_all.MoU$euxinic.MoU[geochem.anox_all.MoU$Mo_U.auth >= eux.cutoff] <- 1
geochem.anox_all.MoU$euxinic.MoU[geochem.anox_all.MoU$Mo_U.auth < eux.cutoff] <- 0

## Bootstrap bimodal euxinia classification using invweight
eux_bootmeans <- as.numeric() 

## Loop through time bins and generate mean authigenic concentration values based on spatial bootstrap
for (bin in (seq(end_age, start_age,  bin_size)+bin_size/2)){
  timebin <- geochem.anox_all.MoU %>%
    filter(time.bin==bin) %>%
    filter(!is.na(euxinic.MoU)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <- one.boot(timebin$euxinic.MoU, mean, R=1000, weights=1/timebin$k)
    eux_bootmeans <- rbind(eux_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}

## Generate dataframe of results and name columns
eux.Mo.U <- as.data.frame(eux_bootmeans)
names(eux.Mo.U) <- c("Prop", "Age")

# Generate an individual bootstrapped mean per time bin
eux.sum.Mo.U <- eux.Mo.U %>%
  group_by(Age) %>%
  summarise(Prop.mean = mean(Prop), Prop.min = mean(Prop) - 2*(sd(Prop, na.rm=T)), Prop.max = mean(Prop) + 2 * (sd(Prop, na.rm=T)))

# Generate panel B plot
eux.plot.MoU <- ggplot(eux.sum.Mo.U, aes(x=Age, y=Prop.mean, ymin=Prop.min, ymax=Prop.max))+
  geom_pointrange(size=1.5)+
  geom_point(shape=21, size=7, stroke=1.5, fill="grey80", alpha=.99)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-0.05,1.05), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression(atop("Proportion ", paste("euxinic (Mo"[auth]*"/U"[auth]*")"))))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.3,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="text", label="B", size=13, y=0.97, x=486)
eux.plot.MoU

# Generate histogram of samples for panel B

geochem.anox_all.MoU$road.riv <- NA
geochem.anox_all.MoU$road.riv[geochem.anox_all.MoU$lithostrat != "Road River Group"] <- "Global Compilation"
geochem.anox_all.MoU$road.riv[geochem.anox_all.MoU$lithostrat == "Road River Group"] <- "Road River Group"

# Summarize counts to ensure histograms span full range of sample numbers
geochem.anox_all.MoU %>% 
  group_by(time.bin) %>%
  tally()

eux.plot.MoU_sample.hist <- ggplot(geochem.anox_all.MoU, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,300), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300), labels=c("","","","150","","","300"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1.5), 
        axis.ticks = element_line(size=1), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position="none",
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
eux.plot.MoU_sample.hist


## ------------------------------------------------------------
# Panel C - proportion samples dominated by authigenic P
## ------------------------------------------------------------

# Authigenic cutoff - choose P concentration. 0.1% or 1000ppm is used as default.
authigenic.cutoff <- 1000

# Filter samples - all anoxic samples with age & P data
geochem.anox_all.P <- filter(geochem, (FeHR.FeT >= 0.38 | FeT.Al >= 0.53 ) & !(is.na(interpreted_age)) & 
                                 !(is.na(P..ppm.)) &
                                 !(is.na(lat_dec)) & !(is.na(long_dec)) &
                                 interpreted_age <= start_age, interpreted_age >= end_age)

## Geospatial weighting coefficients for anoxic shale samples [function "invweight" uses declustering eqn from Keller et al. 2012]
geochem.anox_all.P$k <- invweight(lat = geochem.anox_all.P$lat_dec,
                                    lon = geochem.anox_all.P$long_dec,
                                    age = geochem.anox_all.P$interpreted_age)

# Divide samples into age bins
geochem.anox_all.P$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox_all.P$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2

# Binary coding of authigenic P dominated anoxic samples
geochem.anox_all.P$authigenic.P[geochem.anox_all.P$P..ppm. >= authigenic.cutoff] <- 1
geochem.anox_all.P$authigenic.P[geochem.anox_all.P$P..ppm. < authigenic.cutoff] <- 0

## Bootstrap bimodal authigenic P dominated classification using invweight
P_bootmeans <- as.numeric() 

## Loop through time bins and generate mean proportion values based on spatial bootstrap
for (bin in (seq(end_age, start_age,  bin_size)+bin_size/2)){
    timebin <- geochem.anox_all.P %>%
    filter(time.bin==bin) %>%
    filter(!is.na(authigenic.P)) %>%
    filter(!is.na(k))
  
  if(nrow(timebin)==0){
    print(bin)
  }else{
    bootmean <- one.boot(timebin$authigenic.P, mean, R=1000, weights=1/timebin$k)
    P_bootmeans <- rbind(P_bootmeans, cbind(bootmean$t, rep(bin, nrow(bootmean$t))))
  }}

## Generate dataframe of results and name columns
auth.P <- as.data.frame(P_bootmeans)
names(auth.P) <- c("Prop", "Age")

# Generate an individual bootstrapped mean per time bin
sum.auth.P <- auth.P %>%
  group_by(Age) %>%
  summarise(Prop.mean = mean(Prop), Prop.min = mean(Prop) - 2*(sd(Prop, na.rm=T)), Prop.max = mean(Prop) + 2 * (sd(Prop, na.rm=T)))

# Generate panel C plot
auth.P_plot <- ggplot(sum.auth.P, aes(x=Age, y=Prop.mean, ymin=Prop.min, ymax=Prop.max))+
  geom_pointrange(size=1.5)+
  geom_point(shape=21, size=7, stroke=1.5, fill="grey80", alpha=.99)+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-0.05,0.8), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75))+
  ylab(expression(atop("Proportion ", paste("Authigenic P dominated"))))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.3,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="text", label="C", size=13, y=0.74, x=486)
auth.P_plot

# Generate histogram of samples for anoxic P samples

geochem.anox_all.P$road.riv <- NA
geochem.anox_all.P$road.riv[geochem.anox_all.P$lithostrat != "Road River Group"] <- "Global Compilation"
geochem.anox_all.P$road.riv[geochem.anox_all.P$lithostrat == "Road River Group"] <- "Road River Group"

# Summarize counts to ensure histograms span full range of sample numbers
geochem.anox_all.P %>% 
  group_by(time.bin) %>%
  tally()

anox.plot.P_sample.hist <- ggplot(geochem.anox_all.P, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"))+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,400), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300,350,400), labels=c("","","","","200","","","", "400"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
        axis.line.x = element_blank(), 
        axis.line.y = element_line(lineend = 'square', color="black", size=1.5), 
        axis.ticks = element_line(size=1), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.position="none",
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
anox.plot.P_sample.hist

## ------------------------------------------------------------
# Panel D - Graptolite diversity 
## ------------------------------------------------------------

# Plot includes major evolutionary transitions in the evolution of land plants, based on D'Antonio, Ibarra and Boyce (2020, Geology)

grap.plot <- ggplot(graps, aes(x=bin.center, y=richness))+
  geom_line(size=1)+
  geom_ribbon(aes(ymax=richness, ymin=0), color=NA, fill="grey80")+
  annotate(geom="rect", ymin=110, ymax=118, xmin=-Inf, xmax=470, fill=rgb(166,207,134, maxColorValue = 255), color=NA)+
  annotate(geom="rect", ymin=99, ymax=107, xmin=-Inf, xmax=419.2, fill=rgb(140,176,108, maxColorValue = 255), color=NA)+
  annotate(geom="rect", ymin=88, ymax=96, xmin=-Inf, xmax=393.3, fill=rgb(103,143,102, maxColorValue = 255), color=NA)+
  annotate(geom="text", x=368.5, y=114, size=9, label="Land plants")+
  annotate(geom="text", x=375, y=103, size=9, label="Rooted vascular plants")+
  annotate(geom="text", x=365, y=92, size=9, label="Trees")+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(0,122),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(3, "lines"), unit(3, "lines")),
            size = list(9, 12),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("Graptolite diversity")+xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=36),
        axis.text = element_text( size=32, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="text", label="D", size=13, y=113, x=486)

fig3 <- ggarrange2(eux.plot.Fe_sample.hist, eux.plot.Fe, eux.plot.MoU_sample.hist, eux.plot.MoU, anox.plot.P_sample.hist, auth.P_plot, grap.plot, 
                   ncol=1, heights = c(.2,1,.2,1,.2,1,1))

ggsave(file="Road River Figure 3.pdf", fig3, height=34, width=21)


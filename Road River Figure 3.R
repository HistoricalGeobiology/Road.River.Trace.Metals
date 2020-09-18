library(ggplot2)
library(dplyr)
library(binom)
library(deeptime)
library(egg)

# Set working directory for plot to be saved in
setwd()

# Specify start and end ages
start_age <- 490
end_age <- 360
bin_size <- 10

# Adapt abbreviated stage names in deeptime package for plotting.
stages.renamed <- deeptime::stages
stages.renamed$abbr[69:91] <- c("Fm", "Frs", "Gv", "E", "Em", "P", "Lc", "P", "L", "", "H", "S", "T", "A", "R", "", "Kt", "S", "Drr", "D", "Fl", "Trm", "10")

## ------------------------------------------------------------
# Panel A - proportion euxinic samples based on iron speciation
## ------------------------------------------------------------

# Import data
load("Sperling.Road.River.Data.RData")

# Filter to anoxic samples with ages only
geochem.anox <- filter(geochem, FeHR.FeT >= 0.38 & !(is.na(interpreted_age)))

# Divide samples into age bins
geochem.anox$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2

# Calculate number of samples per bin
n.samples <- as.data.frame(table(geochem.anox$time.bin))

# Calculate number of euxinic samples within each bin
n.successes <-  as.data.frame(table(filter(geochem.anox, Fe.py.FeHR >= 0.7)$time.bin))

# Combine frequency tables of all anoxic and just euxinic samples
binom.input <- cbind(n.samples, n.successes$Freq)

# Rename column headings for clarity
names(binom.input) <- c("Age", "N", "eux")

# Generate dataframe of binomial confidence intervals and means for proportion of euxinic samples
eux.time.sum.Fe <- cbind(binom.input$Age,
                        binom.confint(binom.input$eux, binom.input$N, conf.level  = 0.95, methods = c("asymptotic")))

# Give age column clear name
names(eux.time.sum.Fe)[1] <- "Age"

# Set age as numeric variable
eux.time.sum.Fe$Age <- as.numeric(paste(eux.time.sum.Fe$Age))

# Generate panel A plot
eux.plot.Fe <- ggplot(eux.time.sum.Fe, aes(x=Age, y=mean, ymin=lower, ymax=upper))+
  geom_pointrange(size=1)+
  geom_point(shape=21, size=5, fill="grey80")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-0.05,1.05), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression(atop("Proportion ", paste("euxinic (Fe"[Py]*"/Fe"[HR]*")"))))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.3,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
eux.plot.Fe

# Generate histogram of samples for panel A

geochem.anox$road.riv <- NA
geochem.anox$road.riv[geochem.anox$lithostrat != "Road River Group"] <- "Global Compilation"
geochem.anox$road.riv[geochem.anox$lithostrat == "Road River Group"] <- "Road River Group"

eux.plot.Fe_sample.hist <- ggplot(geochem.anox, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,200), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200), labels=c("","","100", "", "200"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
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
eux.plot.Fe_sample.hist


## ------------------------------------------------------------
# Panel B - proportion euxinic samples based on Mo/U
## ------------------------------------------------------------

# Euxinic cutoff - choose Mo/U ratio. 4 is used as default.
eux.cutoff <- 4

# Filter samples - all anoxic samples with age, Mo, U and Al data
geochem.anox_all.MoU <- filter(geochem, (FeHR.FeT >= 0.38 | FeT.Al >= 0.53 ) & !(is.na(interpreted_age)) & 
                                 !(is.na(Mo..ppm.)) &  !(is.na(U..ppm.)) & !(is.na(Al..wt..)))


# Calculations for authegenic enrichments (Rudnick & Gao for Mo and U)
crustal.Al <- (15.40/(26.981539*2+15.999*3))*(26.981539*2) # calculation from Al2O3 wt% to Al wt%
crustal.Mo <- 1.1
crustal.U <- 2.7
geochem.anox_all.MoU$Mo.det <- (crustal.Mo/crustal.Al) * geochem.anox_all.MoU$Al..wt..
geochem.anox_all.MoU$Mo.auth <- geochem.anox_all.MoU$Mo..ppm. - geochem.anox_all.MoU$Mo.det 
geochem.anox_all.MoU$U.det <- (crustal.U/crustal.Al) * geochem.anox_all.MoU$Al..wt..
geochem.anox_all.MoU$U.auth <- geochem.anox_all.MoU$U..ppm. - geochem.anox_all.MoU$U.det 
geochem.anox_all.MoU$Mo_U.auth <- geochem.anox_all.MoU$Mo.auth/geochem.anox_all.MoU$U.auth 

# Divide samples into age bins
geochem.anox_all.MoU$time.bin <-  seq(end_age, start_age,  bin_size)[as.numeric(cut(geochem.anox_all.MoU$interpreted_age, seq(end_age, start_age,  bin_size)))]+bin_size/2

# Calculate number of samples per bin
n.samples <- as.data.frame(table(geochem.anox_all.MoU$time.bin))

# Calculate number of successes per bin
n.successes <-  as.data.frame(table(filter(geochem.anox_all.MoU, Mo_U.auth >= eux.cutoff)$time.bin))

# Combine frequency tables of all anoxic and just euxinic samples - add zero to end of 485 bin
binom.input <- cbind(n.samples, c(n.successes$Freq, 0))

# Rename column headings for clarity
names(binom.input) <- c("Age", "N", "eux")

# Generate dataframe of binomial confidence intervals and means for proportion of euxinic samples
eux.time.sum.MoU <- cbind(binom.input$Age,
                      binom.confint(binom.input$eux, binom.input$N, conf.level  = 0.95, methods = c("asymptotic")))

# Give age column clear name
names(eux.time.sum.MoU)[1] <- "Age"

# Set age as numeric variable
eux.time.sum.MoU$Age <- as.numeric(paste(eux.time.sum.MoU$Age))

# Generate panel B plot
eux.plot.MoU <- ggplot(eux.time.sum.MoU, aes(x=Age, y=mean, ymin=lower, ymax=upper))+
  geom_pointrange(size=1)+
  geom_point(shape=21, size=5, fill="grey80")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(-0.05,1.05), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab(expression(atop("Proportion ", paste("euxinic (Mo"[auth]*"/U"[auth]*")"))))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(0,1,.3,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(color=NA, fill=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
eux.plot.MoU

# Generate histogram of samples for panel B

geochem.anox_all.MoU$road.riv <- NA
geochem.anox_all.MoU$road.riv[geochem.anox_all.MoU$lithostrat != "Road River Group"] <- "Global Compilation"
geochem.anox_all.MoU$road.riv[geochem.anox_all.MoU$lithostrat == "Road River Group"] <- "Road River Group"

eux.plot.MoU_sample.hist <- ggplot(geochem.anox_all.MoU, aes(time.bin))+
  geom_bar(aes(fill=road.riv))+
  scale_fill_manual(values=c("grey60", "grey30"), name="")+
  coord_cartesian(xlim=c(start_age,end_age-1), ylim=c(0,200), expand=FALSE)+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  scale_y_continuous(breaks=c(0,50,100,150,200), labels=c("","","100", "", "200"))+
  ylab("n")+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,-0.2,1,"cm"),panel.border = element_rect(fill=NA,color=NA, size=2,linetype="solid"),
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
eux.plot.MoU_sample.hist

## ------------------------------------------------------------
# Panel C - Graptolite diversity 
## ------------------------------------------------------------

grap.plot <- ggplot(graps, aes(x=bin.center, y=richness))+
  geom_line(size=1)+
  geom_ribbon(aes(ymax=richness, ymin=0), color=NA, fill="grey80")+
  coord_geo(xlim=c(start_age,end_age-1), ylim=c(0,122),expand=FALSE,
            pos = as.list(rep("bottom", 2)),
            abbrv=list(TRUE, FALSE),
            dat = list(stages.renamed, "periods"),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size = list(5, 8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))+
  theme_minimal()+scale_x_reverse(breaks=c(360,400,440,480))+
  ylab("Graptolite diversity")+xlab("Time (Ma)")+
  theme(plot.margin = margin(.3,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=26),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.position="none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

fig3 <- ggarrange2(eux.plot.Fe_sample.hist, eux.plot.Fe, eux.plot.MoU_sample.hist, eux.plot.MoU, grap.plot, ncol=1, heights = c(.2,1,.2,1,1))

ggsave(file="Road River Figure 3.pdf", fig3, height=20, width=11)


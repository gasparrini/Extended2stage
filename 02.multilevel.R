################################################################################
# R code for the analysis in:
#
#  Sera F, Gasparrini A. Extended two-stage designs for environmental research.
#    Environmental Health. 2022;21:41.
#  https://doi.org/10.1186/s12940-022-00853-z
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/Extended2stage
################################################################################

################################################################################
# MULTILEVEL ANALYSIS BY CITY AND STATE
################################################################################

# LOAD THE PACKAGES
library(mixmeta) ; library(Epi) ; library(maps) ; library(scales)
library(ggplot2)

# LOAD COEF/VCOV FROM FIRST-STAGE MODELS (FOR OZONE)
o3par <- read.csv(file="data/o3par.csv")
coef <- o3par[,grep("coef", names(o3par))]
vcov <- o3par[,grep("vcov", names(o3par))]

# LOAD CITY-LEVEL METADATA AND LINK WITH LAT/LONG DATA
cityinfo <- o3par[,1:4,]
latlong <- read.csv("data/latlong.csv")
cityinfo <- merge(cityinfo, latlong, by="city")

################################################################################
# RUN MODELS

# STANDARD META-ANALYSIS (SINGLE RANDOM-EFFECTS LEVEL)
model0 <- mixmeta(coef, vcov)
summary(model0)

# ADD STATE AS FIXED EFFECTS (PARAMETERIZED TO SHOW STATE-SPECIFIC EFFECTS)
model1 <- mixmeta(coef~0+state, vcov, data=cityinfo)
summary(model1)

# TWO-LEVEL RANDOM EFFECTS (CITY NESTED WITHIN STATE)
model2 <- mixmeta(coef, vcov, data=cityinfo, random=~1|state/city)
summary(model2)

# TEST HETEROGENEITY
qtest(model0)
qtest(model1)
qtest(model2)

# SIGNIFICANCE TEST OF BETWEEN-STATE DIFFERENCES
# NB: REQUIRES STANDARD PARAMETERIZATION AND USE OF ML ESTIMATOR FOR COMPARISON
drop1(update(model1, .~state, method="ml"), test="Chisq")

# ESTIMATED EFFECTS IN RR SCALE FOR 10-UNIT INCREASE
# - POOLED EFFECTS FOR MODELS WITH NO FIXED EFFECTS (WITH CI)
# - RANGE OF EFFECTS FOR MODEL WITH FIXED EFFECTS
ci.exp(model0)
ci.exp(model2)
range(exp(predict(model1)))

################################################################################
# PREDICTIONS

# FIXED-EFFECTS PREDICTION OF STATE-AVERAGE EFFECTS 
statefix <- exp(unique(predict(model1, ci=T)))

# BLUPS AT CITY AND STATE LEVEL FROM TWO-LEVEL MODEL
cityblup <- exp(blup(model2, pi=T))
stateblup <- exp(unique(blup(model2, pi=T, level=1)))

# NAMES
rownames(cityblup) <- cityinfo$cityname[!is.na(coef)]
rownames(stateblup) <- rownames(statefix) <- 
  unique(cityinfo$statename[!is.na(coef)])

################################################################################
# MAP OF CITY-SPECIFIC BLUPS

# EFFECTS, COLOURS, COORDINATES (NB: REVERSE LONG)
cutoff <- pretty(cityblup[,1], 8)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
rrcat <- cut(cityblup[,1], cutoff,  labels=labmap)
colmap <- colorRampPalette(c("yellow","red"))(length(cutoff))
lat <- as.numeric(as.character(cityinfo$lat[!is.na(coef)]))
long <- -as.numeric(as.character(cityinfo$long[!is.na(coef)]))

# MAP OF CITY-SPECIFIC BLUPS (VERSION 1)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=1.8)
legend("bottomright", labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
  title="RR", inset=0.02)
title("Map of risk associated to ozone")
mtext("From a multi-level meta-analysis")

# MAP OF CITY-SPECIFIC BLUPS (VERSION 2)
mapstate <- map_data("state")
ggplot(mapstate, aes(long, lat, group=group)) +
  geom_polygon(fill=NA, col="grey") + 
  borders("usa", col="black") +
  geom_point(aes(group=NULL, col=rrcat), alpha=0.6, size=5, 
    data=data.frame(long=long,lat=lat)) +
  scale_color_brewer(name="RR", palette="YlOrRd") +
  xlim(min(mapstate$long), max(mapstate$long)) +
  ylim(min(mapstate$lat), max(mapstate$lat)) +
  theme_void() +
  theme(legend.position=c(1,0), legend.justification=c(1.2,-0.1)) +
  labs(title="Map of risk associated to ozone", 
     subtitle="From a multi-level meta-analysis") +
  coord_quickmap()

################################################################################
# FOREST PLOT

# SEQUENCES AND LABELS
yseq <- seq(nrow(stateblup))
ylab <- rownames(stateblup)

# FOREST PLOT (VERSION 1)
par(mar=c(5,8.2,4,2)+0.1)
plot(yseq, type="n", xlim=c(0.98,1.02), yaxt="n", xlab="RR", ylab="",
  mgp=c(2.5,1,0), main="Comparison of state-specific risk estimates")
axis(2, at=yseq, las=1, tick=F, labels=ylab, cex.axis=0.8)
grid()
abline(v=1)
arrows(stateblup[,2], yseq-0.2, stateblup[,3], yseq-0.2, col=alpha(2,0.5),
  code=3, angle=90, length=0.01, lwd=2)
points(stateblup[,1], yseq-0.2, col=2, pch=19, cex=1.2)
arrows(statefix[,2], yseq+0.2, statefix[,3], yseq+0.2, col=alpha(4,0.5),
  code=3, angle=90, length=0.01, lwd=2)
points(statefix[,1], yseq+0.2, col=4, pch=19, cex=1.2)
legend("topright", c("Fixed-effects","BLUPs"), pch=19, col=c(4,2), pt.cex=1.2,
  bty="n")
par(mar=c(5,4,4,1)+0.1)

# CREATE DATAFRAME
eststate <- data.frame(rbind(statefix, stateblup), row.names=NULL)
eststate$state <- rep(rownames(statefix), 2)
eststate$model <- rep(c("Fixed-effects","BLUPs"), each=nrow(statefix))
eststate$model <- factor(eststate$model, levels=unique(eststate$model))

# FOREST PLOT (VERSION 2)
ggplot(eststate, aes(fit, state, col=model)) + 
  geom_vline(xintercept=1) +
  geom_errorbar(aes(xmin=ci.lb, xmax=ci.ub), width=0.6, alpha=0.5,
    position=position_dodge(width=0.75)) +
  geom_point(size=2, position=position_dodge(width=0.75)) +
  scale_color_manual(values=c(4,2), name="Prediction")+
  theme_bw() +
  scale_y_discrete(limits=unique(eststate$state)) + 
  coord_cartesian(xlim=c(0.98,1.02)) + 
  theme(legend.position="top") +
  labs(title="Comparison of state-specific risk estimates") +
  xlab("RR") + ylab("")

################################################################################
# SAVE ARTICLE-STYLE PLOTS

# MAP OF CITY-SPECIFIC BLUPS (VERSION 2)
ggplot(mapstate, aes(long, lat, group=group)) +
  geom_polygon(fill=NA, col="grey") + 
  borders("usa", col="black") +
  geom_point(aes(group=NULL, col=rrcat), alpha=0.6, size=5, 
    data=data.frame(long=long,lat=lat)) +
  scale_color_brewer(name="RR", palette="YlOrRd") +
  xlim(min(mapstate$long), max(mapstate$long)) +
  ylim(min(mapstate$lat), max(mapstate$lat)) +
  theme_void() +
  theme(legend.position=c(1,0), legend.justification=c(1.2,-0.1)) +
  coord_quickmap()
ggsave("figures/multilevmap.pdf", height=5, width=9)

# FOREST PLOT (VERSION 2)
ggplot(eststate, aes(fit, state, col=model)) + 
  geom_vline(xintercept=1) +
  geom_errorbar(aes(xmin=ci.lb, xmax=ci.ub), width=0.6, alpha=0.5,
    position=position_dodge(width=0.75)) +
  geom_point(size=2, position=position_dodge(width=0.75)) +
  scale_color_manual(values=c(4,2), name="Prediction")+
  theme_bw() +
  scale_y_discrete(limits=unique(eststate$state)) + 
  coord_cartesian(xlim=c(0.98,1.02)) + 
  theme(legend.position="top") +
  xlab("RR") + ylab("")
ggsave("figures/multilevforest.pdf", height=9, width=5)

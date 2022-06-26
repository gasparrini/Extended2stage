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
# LONGITUDINAL ANALYSIS OF EFFECT MODIFICATION BY AIR CONDITIONING
################################################################################

# LOAD PACKAGES
library(mixmeta) ; library(dlnm) ; library(splines)
library(nlme) ; library(scales)

# LOAD FIRST-STAGE DATA
tmeanperpar <- read.csv(file="data/tmeanperpar.csv")

# LOAD AIR CONDITIONING DATA, RESHAPE TO LONG, REMOVE MISSING
acdata <- read.csv(file="data/acdata.csv")
acdata <- reshape(acdata, varying=seq(ncol(acdata))[-1], idvar="city", sep="",
  timevar="year", direction="long")
acdata <- na.omit(acdata[with(acdata, order(city, year)),])

# SUBSET ORIGINAL DATA TO CITIES WITH AC MEASURES
tmeanperpar <- tmeanperpar[tmeanperpar$city %in% acdata$city,]

# ESTRACT COEF/VCOV FROM FIRST-STAGE MODELS
coef <- as.matrix(tmeanperpar[,grep("coef", names(tmeanperpar))])
vcov <- as.matrix(tmeanperpar[,grep("vcov", names(tmeanperpar))])

# CITY-SPECIFIC META-DATA
cityinfo <- tmeanperpar[,1:6,]

# PERFORM RANDOM-EFFECTS MODEL TO OBTAIN SMOOTH FIT OF AC TREND
mlme <- lme(ac ~ bs(year, degree=2, df=4), data=acdata,
  random=list(city=pdDiag(~ bs(year, degree=2, df=4))))

# PREDICT AC FOR EACH PERIOD IN EACH CITY
cityinfo$ac <- predict(mlme, cityinfo)
summary(cityinfo$ac)
# PROBLEM: ABOVE 100%

################################################################################

# MODEL WITH NO META-PREDICTOR
model0 <- mixmeta(coef, vcov, data=cityinfo, method="ml",
  random=~1|city, bscov="diag")

# MODEL WITH AC AND TIME (TAKES SOME TIME)
# NB: ADD control=list(showiter=TRUE) TO INSPECT THE ITERATIVE OPTIMIZATION
model1 <- mixmeta(coef ~ac+ns(year, knots=1995), vcov, data=cityinfo,
  method="ml", random=~ns(year, knots=1995)|city, bscov="diag")

# TEST (TAKES EVEN LONGER)
drop1(model1, test="Chisq")

################################################################################
# PLOT THE AVERAGE EXPOSURE-RESPONSE AT LOW HIGH AC IN YEAR 2000

# LOAD AVERAGE TEMPERATURE DISTRIBUTION ACROSS CITIES
avgtmeansum <- read.csv("data/avgtmeansum.csv")
tmean <- avgtmeansum$tmean

# DEFINE SPLINE TRANSFORMATION ORIGINALLY USED IN FIRST-STAGE MODELS
knots <- tmean[avgtmeansum$perc %in% paste0(c(50,90), ".0%")]
bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)

# DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
cen <- tmean[which.min(bvar%*%coef(model0))]

# PLOTTING LABELS
xperc <- c(0,1,5,25,50,75,90,99,100)
xval <- tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]

# DEFINE THE VALUES (AC AT 20-80%, YEAR 2000)
datapred <- data.frame(ac=c("Low AC use"=20,"High AC use"=80), year=2000)

# PREDICT COEF/VCOV
pred <- predict(model1, datapred, vcov=T)

# PREDICT ASSOCIATIONS
cp1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
  model.link="log", at=tmean, cen=cen)
cp2 <- crosspred(bvar,coef=pred[[2]]$fit,vcov=pred[[2]]$vcov,
  model.link="log", at=tmean, cen=cen)

# EXPOSURE-RESPONSE PLOT
plot(cp1, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), cex.axis=0.8, col=3,
  ci.arg=list(col=alpha(3,0.3)), main="Temperature-mortality association")
axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=0.8)
lines(cp2, lwd=1.5, col=4, ci="area", ci.arg=list(col=alpha(4,0.3)))
abline(v=cen, lty=2, col=grey(0.8))
mtext("By AC prevalence")
legend("topleft", c(rownames(datapred)), lwd=1.5, col=c(3,4), bty="n", inset=0.1)

################################################################################
# PLOT THE TREND IN RR AT 99TH TMEAN PERC UNDER DIFFERENT AC SCENARIOS

# PREDICT TREND IN AVERAGE AC USE (NB: USE level ARGUMENT IN predict)
acpred <- predict(mlme, data.frame(year=1987:2000), level=0)

# AVERAGE TMEAN AT 99TH PERCENTILE
tmean99 <- tmean[avgtmeansum$perc=="99.0%"]

# DEFINE THE SCENARIOS: CONSTANT AC IN 1987 AND ACTUAL AC TREND
datapred1 <- data.frame(ac=acpred[1], year=1987:2000)
datapred2 <- data.frame(ac=acpred, year=1987:2000)

# PREDICT COEF/VCOV
pred1 <- predict(model1, datapred1, vcov=T)
pred2 <- predict(model1, datapred2, vcov=T)

# PREDICT ASSOCIATIONS
rr99 <- t(sapply(seq(nrow(datapred1)), function(i) {
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp1 <- crosspred(bvar, coef=pred1[[i]]$fit, vcov=pred1[[i]]$vcov,
    model.link="log", at=tmean99, cen=cen)
  cp2 <- crosspred(bvar,coef=pred2[[i]]$fit,vcov=pred2[[i]]$vcov,
    model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp1, c(allRRfit, allRRlow, allRRhigh)), 
    with(cp2, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c(t(outer(c("rr1","rr2"), c("","low","high"), paste0)))
  
  # RETURN
  est
}))

# PLOT
plot(1987:2000, seq(1987:2000), type="n", ylim=range(rr99)*c(0.93,1.07), 
  ylab="RR", xlab="Year", bty="l", las=1, mgp=c(2.5,1,0), cex.axis=0.8,
  main="Trend in risk")
arrows(1987:2000-0.1, rr99[,2], 1987:2000-0.1, rr99[,3], col=alpha(4,0.5),
  code=3, angle=90, length=0.05, lwd=2)
points(1987:2000-0.1, rr99[,1],  type="o", col=4, pch=19)
arrows(1987:2000+0.1, rr99[,5], 1987:2000+0.1, rr99[,6], col=alpha(2,0.5),
  code=3, angle=90, length=0.05, lwd=2)
points(1987:2000+0.1, rr99[,4], type="o", col=2, pch=19)
abline(h=1)
mtext("By scenario of trends in AC prevalence")
legend("top", c("Constant at 1987","Actual trend"), pch=19, col=c(4,2), bty="n",
  inset=0.1)

################################################################################
# SAVE ARTICLE-STYLE PLOT

# GRAPHICALS PARAMETERS
layout(t(1:2))
oldpar <- par(no.readonly = TRUE)
par(mar=c(4,4,2,0.5), cex.axis=0.8)

# PLOTS
plot(cp1, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), col=3,
  ci.arg=list(col=alpha(3,0.3)), main="Temperature-mortality association")
axis(1, at=xval, labels=paste0(xperc, "%"))
lines(cp2, lwd=1.5, col=4, ci="area", ci.arg=list(col=alpha(4,0.3)))
abline(v=cen, lty=2, col=grey(0.8))
legend("topleft", c(rownames(datapred)), lwd=1.5, col=c(3,4), bty="n",
  inset=0.1, title="AC prevalence")

plot(1987:2000, seq(1987:2000), type="n", ylim=range(rr99)*c(0.93,1.07), 
  ylab="RR", xlab="Year", bty="l", las=1, mgp=c(2.5,1,0), main="Trend in risk")
arrows(1987:2000-0.1, rr99[,2], 1987:2000-0.1, rr99[,3], col=alpha(4,0.5),
  code=3, angle=90, length=0.05, lwd=2)
points(1987:2000-0.1, rr99[,1],  type="o", col=4, pch=19)
arrows(1987:2000+0.1, rr99[,5], 1987:2000+0.1, rr99[,6], col=alpha(2,0.5),
  code=3, angle=90, length=0.05, lwd=2)
points(1987:2000+0.1, rr99[,4], type="o", col=2, pch=19)
abline(h=1)
legend("top", c("Constant at 1987","Actual trend"), pch=19, col=c(4,2), bty="n",
  inset=0.1, title="Scenario of AC prevalence")

# RESET
par(oldpar)
layout(1)

# PRINT
#dev.print(pdf, file="figures/longitudinal.pdf", height=5, width=12)

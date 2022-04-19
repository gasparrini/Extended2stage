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
# MODELLING REPEATED MEASURES AND DOSE-RESPONSE RELATIONSHIPS
################################################################################

# LOAD PACKAGES
library(mixmeta) ; library(dlnm) ; library(splines) ; library(scales)

# LOAD COEF/VCOV FROM FIRST-STAGE MODELS
tmeanagepar <- read.csv(file="data/tmeanagepar.csv")
coef <- as.matrix(tmeanagepar[,grep("coef", names(tmeanagepar))])
vcov <- as.matrix(tmeanagepar[,grep("vcov", names(tmeanagepar))])

# CITY-SPECIFIC META-DATA
cityinfo <- tmeanagepar[,1:5,]

################################################################################
# RUN MODELS

# MODEL IGNORING CLUSTERING
model0 <- mixmeta(coef~agecat, vcov, data=cityinfo, method="ml")
summary(model0)

# MODEL CONSIDERING CLUSTERING
model1 <- update(model0, random=~1|city)
summary(model1)

# DEFINE AGE AS CONTINUOUS VARIABLE
cityinfo$age <- c(60,70,85)[factor(cityinfo$agecat, c("under65","65to74","75p"))]

# MODEL WITH LINEAR EFFECT OF AGE
model2 <- update(model1, .~age)
summary(model2)

# MODEL WITH NON-LINEAR EFFECT OF AGE
model3 <- update(model1, .~ns(age, df=2))
summary(model3)

# MODEL COMPARISON AND TESTS
AIC(model0, model1, model2, model3)
BIC(model0, model1, model2, model3)
drop1(model0, test="Chisq")
drop1(model1, test="Chisq")
drop1(model2, test="Chisq")
drop1(model3, test="Chisq")

################################################################################
# PREDICT AND PLOT FOR GIVEN AGES  

# DEFINE THE AGE VALUES
datapred <- data.frame(age=6:9*10)

# PREDICT COEF/VCOV
pred <- predict(model3, datapred, vcov=T)

# LOAD AVERAGE TEMPERATURE DISTRIBUTION ACROSS CITIES
avgtmeansum <- read.csv("data/avgtmeansum.csv")
tmean <- avgtmeansum$tmean

# DEFINE SPLINE TRANSFORMATION ORIGINALLY USED IN FIRST-STAGE MODELS
knots <- tmean[avgtmeansum$perc %in% paste0(c(50,90), ".0%")]
bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)

# PREDICT FOR 80-YEAR-OLD
cen <- tmean[which.min((bvar%*%pred[[3]]$fit))]
cp <- crosspred(bvar, coef=pred[[3]]$fit, vcov=pred[[3]]$vcov,
  model.link="log", at=tmean, cen=cen)

# PLOTTING LABELS
xperc <- c(0,1,5,25,50,75,90,99,100)
xval <- tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]

# PLOT
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), cex.axis=0.8, col=3,
  ci="n", main="Temperature-mortality association")
axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=0.8)
abline(v=cen, lty=2, col=grey(0.8))

# ADD PREDICTIONS FOR OTHER AGES
for(i in c(1,2,4)) {
  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
    model.link="log", at=tmean, cen=cen)
  lines(cp, lwd=1.5, col=i)
}
legend("top", as.character(datapred$age), lwd=1.5, col=1:4, bty="n", inset=0.1,
  title="Age", ncol=4)
mtext("Age patterns")

################################################################################
# SAVE ARTICLE-STYLE PLOT

# GRAPHICALS PARAMETERS
oldpar <- par(no.readonly = TRUE)
par(mar=c(4,4,1,0.5), cex.axis=0.8)

# PREDICT FOR 80-YEAR-OLD
cp <- crosspred(bvar, coef=pred[[3]]$fit, vcov=pred[[3]]$vcov,
  model.link="log", at=tmean, cen=cen)

# PLOT
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), col=3, ci="n")
axis(1, at=xval, labels=paste0(xperc, "%"))
abline(v=cen, lty=2, col=grey(0.8))

# ADD PREDICTIONS FOR OTHER AGES
for(i in c(1,2,4)) {
  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
    model.link="log", at=tmean, cen=cen)
  lines(cp, lwd=1.5, col=i)
}
legend("top", as.character(datapred$age), lwd=1.5, col=1:4, bty="n", inset=0.1,
  title="Age", ncol=4)

# RESET
par(oldpar)

# PRINT
dev.print(pdf, file="figures/dosresp.pdf", height=4.5, width=6)

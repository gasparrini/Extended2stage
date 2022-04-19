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
# POOLING COMPLEX MULTI-PARAMETER ASSOCIATIONS
################################################################################

# LOAD PACKAGES
library(mixmeta) ; library(dlnm) ; library(scales)

# LOAD COEF/VCOV FROM FIRST-STAGE MODELS
tmeanpar <- read.csv(file="data/tmeanpar.csv")
coef <- as.matrix(tmeanpar[,grep("coef", names(tmeanpar))])
vcov <- as.matrix(tmeanpar[,grep("vcov", names(tmeanpar))])

# LINK WITH CENSUS DATA
cityind <- tmeanpar[,1:4,]
citycensus <- read.csv("data/citycensus.csv")
cityind <- merge(cityind, citycensus, by="city")

################################################################################
# RUN THE MODELS

# MODEL WITH NO META-PREDICTOR
model0 <- mixmeta(coef~1, vcov, data=cityind, method="ml")

# SUMMARY AND HETEROGENEITY TEST
summary(model0)
qtest(model0)

# MODELS WITH A SINGLE META-PREDICTOR
# PREDICTORS: TOTAL POP, % OF PEOPLE WITH HIGH-SCHOOL DEGREE, % OF UNEMPLOYED
model1 <- update(model0, .~pop100)
model2 <- update(model0, .~Phigh)
model3 <- update(model0, .~Punem)

# FULL MODEL
model4 <- update(model0, .~pop100+Phigh+Punem)
summary(model4)

# MODEL COMPARISON AND TESTS
AIC(model0, model1, model2, model3, model4)
BIC(model0, model1, model2, model3, model4)
drop1(model4, test="Chisq")

# MODEL SELECTION (STEP FORWARD)
step(model0, .~pop100+Phigh+Punem)

################################################################################
# PLOT THE AVERAGE EXPOSURE-RESPONSE RELATIONSHIPS

# LOAD AVERAGE TEMPERATURE DISTRIBUTION ACROSS CITIES
avgtmeansum <- read.csv("data/avgtmeansum.csv")
tmean <- avgtmeansum$tmean

# DEFINE SPLINE TRANSFORMATION ORIGINALLY USED IN FIRST-STAGE MODELS
knots <- tmean[avgtmeansum$perc %in% paste0(c(50,90), ".0%")]
bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)

# DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
cen <- tmean[which.min(bvar%*%coef(model0))]

# PREDICT THE ASSOCIATION
cp <- crosspred(bvar, coef=coef(model0), vcov=vcov(model0), model.link="log",
  at=tmean, cen=cen)

# PLOTTING LABELS
xperc <- c(0,1,5,25,50,75,95,99,100)
xval <- tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]

# PLOT
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, xaxt="n", mgp=c(2.5,1,0), cex.axis=0.8,
  main="Temperature-mortality association")
axis(1, at=xval, labels=paste0(xperc, "%"))
abline(v=cen, lty=2, col=grey(0.8))

################################################################################
# PREDICT AND PLOT FOR GIVEN VALUES OF META-PREDICTORS  

# DEFINE THE VALUES (SMALL/LARGE POP, WITH AVERAGE OF OTHERS)
datapred <- data.frame(
  pop100 = quantile(cityind$pop100, c(5,95)/100),
  Phigh = mean(cityind$Phigh),
  Punem = mean(cityind$Punem),
  row.names=c("Small population", "Large population")
)

# PREDICT COEF/VCOV
pred <- predict(model4, datapred, vcov=T)

# PREDICT ASSOCIATIONS
cp1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
  model.link="log", at=tmean, cen=cen)
cp2 <- crosspred(bvar,coef=pred[[2]]$fit,vcov=pred[[2]]$vcov,
  model.link="log", at=tmean, cen=cen)

# PLOT
plot(cp1, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), cex.axis=0.8, col=3,
  ci.arg=list(col=alpha(3,0.3)), main="Temperature-mortality association")
axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=0.8)
lines(cp2, lwd=1.5, col=4, ci="area", ci.arg=list(col=alpha(4,0.3)))
abline(v=cen, lty=2, col=grey(0.8))
mtext("By population size")
legend("topleft", c(rownames(datapred)), lwd=1.5, col=c(3,4), bty="n", inset=0.1)

################################################################################
# SAVE ARTICLE-STYLE PLOT

# GRAPHICALS PARAMETERS
layout(t(1:2))
oldpar <- par(no.readonly = TRUE)
par(mar=c(4,4,2,0.5), cex.axis=0.8)

# PLOTS
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, xaxt="n", mgp=c(2.5,1,0))
axis(1, at=xval, labels=paste0(xperc, "%"))
abline(v=cen, lty=2, col=grey(0.8))
title("Average relationship")

plot(cp1, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
  lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), col=3,
  ci.arg=list(col=alpha(3,0.3)))
axis(1, at=xval, labels=paste0(xperc, "%"))
lines(cp2, lwd=1.5, col=4, ci="area", ci.arg=list(col=alpha(4,0.3)))
abline(v=cen, lty=2, col=grey(0.8))
legend("topleft", c(rownames(datapred)), lwd=1.5, col=c(3,4), bty="n", inset=0.1)
title("By population size")

# RESET
par(oldpar)
layout(1)

# PRINT
dev.print(pdf, file="figures/multipar.pdf", height=5, width=12)

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
# RUN FIRST-STAGE MODELS AND SAVE DATA
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(mixmeta)
library(tsModel) ; library(splines) 
library(lubridate)

# PATH TO NMMAPS CITY-SPECIFIC SERIES
# NB: TO BE REPLACED BY THE PATH TO THE ORIGINAL NMMAPS DATA
path <- ""

# LOAD THE DATASET WITH METADATA FOR THE CITIES
cities <- read.csv("data/cities.csv")

################################################################################
# LOOP

# CREATE OBJECTS TO STORE THE RESULTS
tmeanparlist <- tmeanperparlist <- tmeanageparlist <- o3parlist <-
  tmeansumlist <- vector("list", nrow(cities))

# RUN THE LOOP
for(i in seq(nrow(cities))) {
  
  # PRINT CITY
  cat(cities$city[i],"")
  
  # LOAD THE DATA NOT COLLAPSED BY AGE
  load(paste0(path, "/cities/", cities$city[i]))
  
  # SELECT VARIABLES
  vars <- c("date","dow","death","accident","tmpd","o3tmean","o3mtrend","agecat")
  data <- subset(data, select=vars)
  
  # CONVERT ENVIRONMENTAL VARIABLES
  # - TEMPERATURE FROM FAHRENHEIT TO CELSIUS
  # - OZONE FROM DETRENDED AND CONVERTED IN MICROGR/M3
  data <- transform(data,
    tmean = (tmpd-32)*5/9,
  	o3 = (o3tmean+o3mtrend) * 1.96
  )

  # CREATE COLLAPSED DATA (SELECT ONE AGE GROUP AND COLLAPSE DEATH CAUSES)
  dataggr <- data[data$agecat=="under65",]
  dataggr$death <- tapply(data$death, data$date, sum)
  dataggr$accident <- tapply(data$accident, data$date, sum)
  
  # CREATE ALL-CAUSE
  data$all <- data$death + data$accident
  dataggr$all <- dataggr$death + dataggr$accident
  
  # SUBSET FOR SUMMER-ONLY
  datasum <- subset(data, month(date) %in% 6:9)
  dataggrsum <- subset(dataggr, month(date) %in% 6:9)

################################################################################
# ANALYSIS OF TEMPERATURE - ALL-CAUSE MORTALITY (SUMMER-ONLY)
	
	# DEFINE CROSS-BASIS FOR TEMPERATURE
  cbtmean <- crossbasis(dataggrsum$tmean, lag=3, argvar=list(fun="bs", degree=2,
    knots=quantile(dataggrsum$tmean, c(50,90)/100, na.rm=T)), 
    arglag=list(fun="integer"), group=year(dataggrsum$date))
  
  # RUN THE MODEL AND EXTRACT REDUCED PRED (RE-CENTRED LATER)
  model <- glm(all ~ cbtmean + dow + ns(yday(date), df=4)*factor(year(date)),
    data=dataggrsum, family=quasipoisson)
  redpred <- crossreduce(cbtmean, model, cen=mean(dataggrsum$tmean, na.rm=T))
  
  # STORE PARAMETERS (COEF + VECTORIZED VCOV)
  ncoef <- length(coef(redpred))
  par <- c(coef(redpred), vechMat(vcov(redpred)))
  names(par) <- c(paste0("coef", seq(ncoef)),
    paste0("vcov", seq(ncoef*(ncoef+1)/2)))
  tmeanpar <- data.frame(cities[i, c("city", "cityname", "state", "statename")],
    t(par), row.names=i)
  tmeanparlist[[i]] <- tmeanpar

################################################################################
# ANALYSIS OF TEMPERATURE - ALL-CAUSE MORTALITY (SUMMER-ONLY, BY PERIOD)
  
  # DEFINE THE PERIODS
  yearlist <- list(1987:1989, 1990:1992, 1993:1995, 1996:1998, 1999:2000)
  
  # PERFORM MODEL BY PERIOD
  parlist <- lapply(yearlist, function(ysub) {
    model <- glm(all ~ cbtmean + dow + ns(yday(date), df=4)*factor(year(date)),
      data=dataggrsum, family=quasipoisson, subset=year(date) %in% ysub)
    redpred <- crossreduce(cbtmean, model, cen=mean(dataggrsum$tmean, na.rm=T))
    t(c(coef(redpred), vechMat(vcov(redpred))))
  })
  
  # STORE PARAMETERS (COEF + VECTORIZED VCOV)
  par <- do.call(rbind, parlist)
  colnames(par) <- names(tmeanpar)[-c(1:4)]
  tmeanperpar <- data.frame(
    cities[i, c("city", "cityname", "state", "statename")],
    period = sapply(yearlist, function(x) paste(range(x), collapse="-")),
    year = sapply(yearlist, mean),
    par,
    row.names=paste0(i, ".", seq(yearlist))
  )
  tmeanperparlist[[i]] <- tmeanperpar
  
################################################################################
# ANALYSIS OF TEMPERATURE - ALL-CAUSE MORTALITY (SUMMER-ONLY, BY AGE)

  # DEFINE AGE GROUPS
  agecat <- as.character(unique(data$agecat))
  
  # PERFORM MODEL BY AGE (SELECTING FROM NON-COLLAPSED DATA)
  parlist <- lapply(agecat, function(cat) {
    y <- data$all[data$agecat==cat & month(data$date) %in% 6:9]
    model <- glm(y ~ cbtmean + dow + ns(yday(date), df=4)*factor(year(date)),
      data=dataggrsum, family=quasipoisson)
    redpred <- crossreduce(cbtmean, model, cen=mean(dataggrsum$tmean, na.rm=T))
    t(c(coef(redpred), vechMat(vcov(redpred))))
  })

  # STORE PARAMETERS (COEF + VCOV FOR 10-UNIT INCREASE)
  par <- do.call(rbind, parlist)
  colnames(par) <- names(tmeanpar)[-c(1:4)]
  tmeanagepar <- data.frame(
    cities[i, c("city", "cityname", "state", "statename")],
    agecat=agecat,
    par,
    row.names=paste0(i, ".", seq(agecat))
  )
  tmeanageparlist[[i]] <- tmeanagepar
  
################################################################################
# ANALYSIS OF OZONE - NON-EXTERNAL MORTALITY (FULL YEAR)
  
  # DEFINE MOVING AVERAGE OF OZONE AT LAG 0-1
  o301 <- runMean(dataggr$o3, 0:1)
  
  # DEFINE CROSS-BASIS FOR TEMPERATURE
  cbtmean <- crossbasis(dataggr$tmean, lag=3, argvar=list(fun="bs", degree=2,
    knots=quantile(dataggrsum$tmean, c(10,50,90)/100, na.rm=T)), 
    arglag=list(fun="strata"))
  
  # RUN THE MODEL AND EXTRACT PAR (ONLY IF ENOUGH NON-MISSING)
  par <- if(nrow(na.omit(cbind(o301, cbtmean))) < 500 ) c(NA,NA) else {
    model <- glm(death ~ o301 + cbtmean + dow + ns(date, df=14*7),
      data=dataggr, family=quasipoisson)
    c(coef(model)["o301"]*10, vcov(model)["o301","o301"]*10)
  }
  names(par) <- c("coef", "vcov")

  # STORE PARAMETERS (COEF + VCOV FOR 10-UNIT INCREASE)
  o3par <- data.frame(cities[i, c("city", "cityname", "state", "statename")],
    t(par), row.names=i)
  o3parlist[[i]] <- o3par
  
################################################################################
# TEMPERATURE DISTRIBUTION (SUMMER ONLY)

  # DEFINE PERCENTILES
  per <- c(0:9/10, 1:99, 991:1000/10)/100
  tmeansumlist[[i]] <- quantile(dataggrsum$tmean, per, na.rm=T)
}

################################################################################
# PREPARE AND STORE

# RBIND COEF/VCOV TOGETHER IN DATAFRAMES
tmeanpar <- do.call(rbind, tmeanparlist)
tmeanperpar <- do.call(rbind, tmeanperparlist)
tmeanagepar <- do.call(rbind, tmeanageparlist)
o3par <- do.call(rbind, o3parlist)

# CREATE COUNTRY-AVERAGE SUMMER TEMPERATURE DISTRIBUTION
avgtmeansum <- data.frame(perc=names(tmeansumlist[[1]]), 
  tmean=apply(do.call(cbind, tmeansumlist), 1, mean))

# STORE
write.csv(tmeanpar, file="data/tmeanpar.csv", row.names=F)
write.csv(tmeanperpar, file="data/tmeanperpar.csv", row.names=F)
write.csv(tmeanagepar, file="data/tmeanagepar.csv", row.names=F)
write.csv(o3par, file="data/o3par.csv", row.names=F)
write.csv(avgtmeansum, file="data/avgtmeansum.csv", row.names=F)

# Extended two-stage designs for environmental research

An illustration of a unified framework that enables several design extension of standard two-stage analyses applied in environmental epidemiology 


--------------------------------------------------------------------------------

This repository stores the updated R code and data to reproduce the four case studies presented in the article:

Sera F, Gasparrini A. Extended two-stage designs for environmental research. *Environmental Health*. 2022;21:41. DOI: 10.1186/s12940-022-00853-z. [[freely available here](http://www.ag-myresearch.com/2022_sera_envhealth.html)]


### Data

The various extensions are illustrated using the database collected as part of the National Morbidity, Mortality and Air Pollution Study (NMMAPS) . This database contains, among other information, daily series of mortality counts as well as weather and pollution measurements for the period 1987–2000 in each of 108 cities in the USA. These data are here complemented with information on air conditioning use, collected longitudinally for a subset of cities (see the article for sources and links).

**Note**: The original datasets composing the NMMAPS database were collected on the 15th of May, 2013, through the R package NMMAPSdata. The package is now archived and the mortality series are not available anymore. Therefore, the results of the first-stage models in each of the four case studies are provided separately in the folder *data*, together with city-specific meta-variables originally in the NMMAPS database. The code for fitting the first-stage models and store the data is provided for completeness (see below), although it cannot be reproduced without the original data.

Specifically, the folder *data* includes the following datasets:

  * *cities.csv*, *citycensus.csv*, *counties.csv*, *latlong.csv*, and *region.csv*: original datasets available in the NMMAPS database, with meta-variables providing information on city-speficic characteristic.
  * *acdata.csv*: complementary data with information on air conditioning use, collected longitudinally for a subset of cities.
  * *tmeanpar.csv*: estimates of the association between heat (summer-only) and all-cause mortality in each city. The estimates are represented by the coefficients and lower triangular entries of the (co)variance matrix of the spline function applied to model summer temperature over lag 0--3. See the Supplementary Information in the article for modelling details.
  * *tmeanagepar.csv* and *tmeanperpar.csv*: estimates of the same association between heat (summer-only) and all-cause mortality in each city, but stratified by age group (0–64, 65–74, 65 and older) and period (1987--98, 1990--92, 1993--95, 1996--98, and 1999--2000), respectively. See the Supplementary Information in the article for modelling details.
  * *o3par.csv*: estimates of the association between ozone and non-accidental mortality in each city. The estimates are represented by the coefficient and variance of log-RR for an increase in ozone of 10 microg/m3. See the Supplementary Information in the article for modelling details.


### R code

The four main R scripts perform the analyses in each of the four case studies. Specifically:

  * *01.heterogeneity.R* reproduces the first case study, illustrating the extension of the two-stage design for pooling multi-parameter associations. Specifically, the case study investigates the non-linear and delayed relationship between heat and all-cause mortality during the summer months and the potential role of city-specific characteristics in modifying the risk.
  * *02.multilevel.R* reproduces the second case study, describing  the extension of the two-stage design for the analysis of complex hierarchical structures and geographical clustering. Specifically, the case study demonstrates how to combine estimates of associations between ozone and mortality in multiple cities nested within states.
  * *03.doseresp.R* performs the analysis of the third case study, showing the extension of the two-stage design for sub-groups analysis and dose–response relationships. Specifically, the case study illustrates how the heat-mortality association can be estimated when there are repeated measures from the same city, resulting from multiple first-stage models fitted by different age groups. The case study also shows how to flexibly model the age effect in a dose-response fashion.
  * *04.longitudinal.R* reproduces the fourth case study, illustrating the extension of the two-stage design for longitudinal analysis of estimates collected along time. Specifically, the case study examines the temporal changes in the exposure–response curve between heat and mortality, and assesses the role of air conditioning (AC) in attenuating the risk.
  

### Additional R code

The additional script *00.prepdata.R* in the folder *addcode* folder was used to create the datasets included in the folder *data*. Specifically, the code loads the original NMMAPS data and performs the first-stage time series models to obtain exposure summaries and estimates of short-term risk associations with air pollution and temperature. **Note**: to run this scripts, the user needs to gather the original NMMAPS data, which is not made available here. The script is made available mostly for information about the first-stage modelling, and it is not meant to be run by the users.

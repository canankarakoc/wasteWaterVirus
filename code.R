# Virus prevalence in waste water treatment plants 
# PCR copy numbers vs. incidence 
# Author: Canan Karakoc

setwd("~/GitHub/wasteWaterVirus")
library(tidyverse)
library(readxl)
library(forecast)
library(zoo)
library(tseries)
library(lubridate)
library(stringr)
library(bsts)
library(vars)

# raw data 
filenames <- list.files("data/raw", pattern="*.xlsx", full.names = F)
dir_path <-"~/GitHub/wasteWaterVirus/data/raw/"

# read all data in the folder and save as a list
datalist <- list()
for (i in 1:length(filenames)){
  data <- read_excel(path = paste0(dir_path, filenames[i]), col_names = T)
  datalist[[length(datalist) + 1 ]] <- data
}
names(datalist) <- filenames

# separate date column 
sepfun      <- function(data) separate(data, date, into = c("week", "year"))
newList     <- lapply(datalist, sepfun)

# convert all to numeric 
numfun   <- function(data) mutate_if(data, is.character, as.numeric)
lastList <- lapply(newList, numfun)

# extrapolate missing data when needed 
extrapolate <- function(x) as.numeric(na.approx(zoo(x)))

###############
# Leipzig sars
###############
leip_sars <- lastList[[1]]
leip_sars$copy <- extrapolate(leip_sars$copy)
ts1_copy <- ts(leip_sars[3], frequency = 1)
ts1_inc  <- ts(leip_sars[4], frequency = 1)
  
acf(ts1_copy)
# there is an upward trend, which means the ts is notstationary
# differentiating might solve it 
diff_ts1_copy <- as.numeric(diff(ts1_copy))
acf(diff_ts1_copy) # now it is stationary
adf.test(diff_ts1_copy) 

# incidence
acf(ts1_inc) # incidence seems more autocorrelated
diff_ts1_inc <- as.numeric(diff(ts1_inc))

# partial autocorrelattion 
# The partial autocorrelation is the amount of correlation between a time 
# series and lags of itself that is not explained by a previous lag
pacf(diff_ts1_copy)
pacf(diff_ts1_inc)

# did I chose the lag differences correctly?
ndiffs(x = ts1_copy) # for copy numbers one lag 
ndiffs(x = ts1_inc)  # increment doesn't need transformation
# but the values are too big maybe not a bad idea

# ARIMA models 
# data is not big, and predictions are not so good so I 
# do not choose training/test data atm 
m1 <- auto.arima(diff_ts1_copy, method='ML') # it has only MA part 
m2 <- auto.arima(diff_ts1_inc, method='ML') # this has MA

summary(m1)
# which one predicts the future better?
predict(m1, n.ahead = 5, se.fit=TRUE)
forecast1 <- forecast(object = m1, h=5) 
plot(forecast1)

forecast2 <- forecast(object = m2, h=5) 
plot(forecast2)

# Maybe they predict each others dynamics
# Vector autoregressive model (VAR) 
multiTS1 <- ts(data = leip_sars[3:4]) 

numDiffs  <- ndiffs(multiTS1)
multiTS1_diff <- diff(multiTS1, differences = numDiffs)
plot(multiTS1_diff)

var1 <- VAR(multiTS1_diff, lag.max = 6)
var1$varresult

predict(var1, n.ahead = 5)
forecast_multi <- forecast(object = var1, h = 5) 
plot(forecast_multi)
# predictions got better in multivariate time series
# if you have both information, you can predict the disease prevalence better
# these models are like linear fit if you look at the model summary 
summary(var1)
AIC(var1) #too big though
# for instance copy numbers best explained by it's own past 3 observation and 
# past three observation of incidence data. R2 is low 
# but it means you can predict with 40% confidence 
# on the other hand, incidence can predict itself with it's own 2 past observation
# we can do this now for all the viruses but first look at the cross correlation
# between the time series 

print(ccf(as.numeric(diff_ts1_copy), as.numeric(diff(ts1_inc))))
# cross correlation is above confidence intervals around lag 1 
# not bad,  copy numbers can almost predict incidences 

# I'll do the same for the others 
# not in a loop 

# all datasets 

freiberg_noro <- lastList[[2]]
leipzig_flu   <- lastList[[3]]
leipzig_astro <- lastList[[4]]
leipzig_rota  <- lastList[[5]]
leipzig_noro  <- lastList[[6]]
freiberg_sars <- lastList[[7]]
freiberg_fluA <- lastList[[8]]
freiberg_astro<- lastList[[9]]
freiberg_rota <- lastList[[10]]

####################
# freiberg norovirus
####################
acf(extrapolate(freiberg_noro$copy_gen1)) # not stationary 
acf(extrapolate(freiberg_noro$copy_gen2)) 
acf(extrapolate(freiberg_noro$copy_gen3)) # not stationary
acf(extrapolate(freiberg_noro$incidence)) # this too 

ndiffs(x = freiberg_noro$copy_gen1)
ndiffs(x = freiberg_noro$copy_gen3)
ndiffs(x = freiberg_noro$incidence)

# I think I'll go ahead and do the multivariate model for all 
# we have to summarise the results somehow 
# maybe show R2 for all models 

# one run arima as above to show incidence is better predicted with qpcr data

multiTS2.1 <- ts(data = na.approx(zoo(freiberg_noro[c(6,3)])))
multiTS2.2 <- ts(data = na.approx(zoo(freiberg_noro[c(6,4)])))
multiTS2.3 <- ts(data = na.approx(zoo(freiberg_noro[c(6,5)])))

ndiffs(multiTS2.1);ndiffs(multiTS2.2);ndiffs(multiTS2.3)
multiTS2.1_diff <- diff(multiTS2.1, differences = 1)
multiTS2.2_diff <- diff(multiTS2.2, differences = 1)
multiTS2.3_diff <- diff(multiTS2.3, differences = 1)

var2.1 <- VAR(multiTS2.1_diff, lag.max = 6); summary(var2.1)
var2.2 <- VAR(multiTS2.2_diff, lag.max = 6); summary(var2.2)
var2.3 <- VAR(multiTS2.3_diff, lag.max = 6); summary(var2.3)

forecast_multi2.1 <- forecast(object = var2.1, h = 5) 
forecast_multi2.2 <- forecast(object = var2.2, h = 5) 
forecast_multi2.3 <- forecast(object = var2.3, h = 5) 

plot(forecast_multi2.1);plot(forecast_multi2.2);plot(forecast_multi2.3)

# cross correlations
# discuss how to summarise 
fn_gen1 <- diff(extrapolate(freiberg_noro$copy_gen1))
fn_gen2 <- diff(extrapolate(freiberg_noro$copy_gen2))
fn_gen3 <- diff(extrapolate(freiberg_noro$copy_gen3))
fn_inc  <- diff(extrapolate(freiberg_noro$incidence))

ccf(fn_gen1, fn_inc) #lag 10 is cross correlated? it's odd or genotype 1 can predict incidence
#10 time points earlier, coefficient is too low so I'd ignore is 
ccf(fn_gen2, fn_inc) # around lag4 positive correlation
ccf(fn_gen3, fn_inc) # same there is some correlation but not too high
# but interesting genotype1 and genotype2 are correlated

#############
# leipzig flu
#############
# this is not a good data overall, too many s0s 
multiTS3  <- ts(data = na.approx(zoo(leipzig_flu[c(4,3)])))
ndiffs(multiTS3)
var3 <- VAR(multiTS3, lag.max = 6); summary(var3)
forecast_multi3 <- forecast(object = var3, h = 5) 
plot(forecast_multi3) #as expected, horrible 

# cross correlations
print(ccf(extrapolate(leipzig_flu$incidence), extrapolate(leipzig_flu$copy)))
# nope, don't trust this although lag 1-2 looks correlated 
# there is non-stationarity and differencing woulndt work 

# leipzig astro
multiTS4  <- ts(data = na.approx(zoo(leipzig_astro[c(4,3)])))
ndiffs(multiTS4)
acf(multiTS4)

#model
var4 <- VAR(multiTS4, lag.max = 6); summary(var4)
forecast_multi4 <- forecast(object = var4, h = 5) 
plot(forecast_multi4) # amazing!
summary(var4) # copy predicts itself, incidence is predicted by itself and copy number lag2-4

# cross correlations
print(ccf(extrapolate(leipzig_astro$copy), extrapolate(leipzig_astro$incidence)))
# lag 0-7 pretty correlated 
# but it is negative lag, so incidence predicts the copy numbers 

###############
# leipzig rota
###############
multiTS5  <- ts(data = na.approx(zoo(leipzig_rota[c(4,3)])))
ndiffs(multiTS5)
acf(multiTS5)

#take the diff
multiTS5_diff <- diff(multiTS5, differences = 1)

#model
var5 <- VAR(multiTS5_diff, lag.max = 6)
forecast_multi5 <- forecast(object = var5, h = 5) 
plot(forecast_multi5) 
summary(var5) # WOW!!! Incidence is literally predicted by 1-3 lags of copy numbers

# cross correlations
lr_copy <- diff(extrapolate(leipzig_rota$copy))
lr_inc  <- diff(extrapolate(leipzig_rota$incidence))
print(ccf(lr_copy, lr_inc)) # there is correlation lag 5-1 but one negative, one positive 
#i dont know if it is random 

###############
# leipzig noro
###############
multiTS6.1  <- ts(data = na.approx(zoo(leipzig_noro[c(6,3)])))
multiTS6.2  <- ts(data = na.approx(zoo(leipzig_noro[c(6,4)])))
multiTS6.3  <- ts(data = na.approx(zoo(leipzig_noro[c(6,5)])))

ndiffs(multiTS6.1);ndiffs(multiTS6.2);ndiffs(multiTS6.3)
acf(multiTS6.1);acf(multiTS6.2);acf(multiTS6.3)

# I want to check the acfs separately to be sure 
acf(extrapolate(leipzig_noro[3]));acf(extrapolate(leipzig_noro[4]));acf(extrapolate(leipzig_noro[5]));acf(extrapolate(leipzig_noro[6]))

acf(diff(extrapolate(leipzig_noro[3]))) #ok, good taking the difference normalizes the genotype1

#take the diff
multiTS6.1_diff <- diff(multiTS6.1, differences = 1)
multiTS6.2_diff <- diff(multiTS6.2, differences = 1)
multiTS6.3_diff <- diff(multiTS6.3, differences = 1)

#model
var6.1 <- VAR(multiTS6.1_diff, lag.max = 6); summary(var6.1) # genotype 1 is not promissing 
forecast_multi6.1 <- forecast(object = var6.1, h = 5) 
plot(forecast_multi6.1) 

var6.2 <- VAR(multiTS6.2_diff, lag.max = 6); summary(var6.2) # seriously?
forecast_multi6.2 <- forecast(object = var6.2, h = 5) 
plot(forecast_multi6.2) 

var6.3 <- VAR(multiTS6.3_diff, lag.max = 6); summary(var6.3) # genotype 1 is not promissing 
forecast_multi6.3 <- forecast(object = var6.3, h = 5) 
plot(forecast_multi6.3) # incidence definitely predict itself 

# cross correlations
ln_gene1_copy <- diff(extrapolate(leipzig_noro$copy_gen1))
ln_gene2_copy <- diff(extrapolate(leipzig_noro$copy_gen2))
ln_gene3_copy <- diff(extrapolate(leipzig_noro$copy_gen3))
ln_inc        <- diff(extrapolate(leipzig_noro$incidence))

print(ccf(ln_gene1_copy, ln_inc)) # as expected, genotpe 1 is not corelated
print(ccf(ln_gene2_copy, ln_inc)) 
print(ccf(ln_gene3_copy, ln_inc)) # others some lags but low

# I think it there is some seasonality effect remained
# sometimes every 5-10 lag is correlated or so 
# but no obvious seasonality 

###############
# freiberg_sars
###############

multiTS7  <- ts(data = na.approx(zoo(freiberg_sars[c(4,3)])))
ndiffs(multiTS7)
acf(multiTS7)

#take the diff
multiTS7_diff <- diff(multiTS7, differences = 1)

#model
var7 <- VAR(multiTS7_diff, lag.max = 6); summary(var7) # MEH...
forecast_multi5 <- forecast(object = var5, h = 5) 
plot(forecast_multi5) 

# cross correlations
fs_copy <- diff(extrapolate(freiberg_sars$copy))
fs_inc  <- diff(extrapolate(freiberg_sars$incidence))
print(ccf(fs_copy, fs_inc)) # meh

################
# freiberg fluA
################
# I have 0 expectations
multiTS8  <- ts(data = na.approx(zoo(freiberg_fluA[c(4,3)])))
ndiffs(multiTS8)
acf(multiTS8)

#model
var8 <- VAR(multiTS8, lag.max = 6)
forecast_multi8 <- forecast(object = var8, h = 5) 
plot(forecast_multi8) 
summary(var5) # well... 

# cross correlations
ff_copy <- diff(extrapolate(freiberg_fluA$copy))
ff_inc  <- diff(extrapolate(freiberg_fluA$incidence))
print(ccf(ff_copy, ff_inc)) # there is correlation, yes but too many 0s to talk about it

###############
# freiberg_astro
###############
multiTS9  <- ts(data = na.approx(zoo(freiberg_astro[c(4,3)])))
ndiffs(multiTS9)
acf(multiTS9)

#take the diff
multiTS9_diff <- diff(multiTS9, differences = 1)

#model
var9 <- VAR(multiTS9_diff, lag.max = 6); summary(var9) 
forecast_multi9 <- forecast(object = var9, h = 5) 
plot(forecast_multi9) # lag 2-3 of copy numbers 

# cross correlations
fa_copy <- diff(extrapolate(freiberg_rota$copy))
fa_inc  <- diff(extrapolate(freiberg_rota$incidence))
print(ccf(fa_copy, fa_inc)) 

###############
# freiberg_rota
###############
multiTS10  <- ts(data = na.approx(zoo(freiberg_rota[c(4,3)])))
ndiffs(multiTS10)
acf(multiTS19)

#take the diff
multiTS10_diff <- diff(multiTS10, differences = 1)

#model
var10 <- VAR(multiTS10_diff, lag.max = 6); summary(var10) 
forecast_multi10 <- forecast(object = var10, h = 5) 
plot(forecast_multi10) # Well, incidence peeks comes before copy numbers. 
# if there is no incidence data available copy number can predict the disease prevalence
# my dumb interpretation

# cross correlations
fr_copy <- diff(extrapolate(freiberg_rota$copy))
fr_inc  <- diff(extrapolate(freiberg_rota$incidence))
print(ccf(fr_copy, fr_inc)) # yeah pretty much copy number follows incidence with lag 1 
# 0.474 correlation coefficient


setwd("/Users/suryalamichhane/Desktop/TS.Proj")
library(plyr)
library(ggplot2)
library(forecast)
library(dlm)
library(astsa) # Stoffer's package
df.pokhara = read.table("df.Pokhara.txt", header = TRUE) 
pokhara.rainfall = read.table("PokharaRainfall.txt", header = TRUE)
colnames(pokhara.rainfall) = NULL
pokhara.rainfall = as.vector(t(pokhara.rainfall[,-1]))
pokhara.rainfall = ts(pokhara.rainfall, frequency = 12, start = c(1980,1))
plot(pokhara.rainfall, type = 'o', pch = 20, ylab = ~mm,
main = "Pokhara rainfall 1980-2015  monthly average")
lines(lowess(pokhara.rainfall, f = 0.3), lwd = 2, col = "hotpink") ## add smoothed trend
plot(stl(pokhara.rainfall, s.window = "periodic"))
fit.rain <- StructTS(pokhara.rainfall, "level")
fit.rain
lines(fitted(fit.rain), lty = "dashed", lwd = 2)
lines(tsSmooth(fit.rain), lty = "dotted", lwd = 2)
tsdiag(fit.rain)
####
colnames(df.pokhara) = c("Month", "T_max", "T_min")
pokhara.temp = ts(df.pokhara[, -1], frequency = 12, start = c(1980,1))
pokhara.max =  pokhara.temp[, -2]
pokhara.min = pokhara.temp[, -1]
# maximum and min plot
par(mfrow=c(3, 1))
plot(pokhara.max, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Pokhara air temperature max from 1980-2015 average")
lines(lowess(pokhara.max, f = 0.3), lwd = 2, col = "hotpink") ## add smoothed trend
#Minimum
plot(pokhara.min, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Pokhara air temperature min from 1980-2015 average")
lines(lowess(pokhara.min, f = 0.3), lwd = 2, col = "hotpink") ## add smoothed trend
## check the trend and seasonality of maximum
plot(stl(pokhara.max, s.window = "periodic"))
## check the trend and seasonality of minimun
plot(stl(pokhara.min, s.window = "periodic"))
##
#####
###### Particulate and mortality in LA (in package astsa)
######

## basic
plot(pokhara.temp[, c("T_max", "T_min")], main = "Max and Min Temperature of ")
##
## a bit fancier
T_max = pokhara.temp[,1]
T_min = pokhara.temp[,2]
opar <- par(mfrow=c(2,1), mar=c(0,3.5,0,3), oma=c(3.5,0,2,0), mgp=c(2,.6,0), cex.lab=1.1, tcl=-.3, las=1)
plot(T_max, ylab=expression(T[~max]~~~~(degree~C)), xaxt="no", type='n')
  grid(lty=1, col=gray(.9))
  lines(T_max, col=rgb(0,0,.9))
plot(T_min, ylab=expression(T[~min]~~~~(degree~C)), xaxt="no", yaxt='no', type='n')
  grid(lty=1, col=gray(.9))
  lines(T_min, col=rgb(.9,0,.9)) 
  axis(4) 
title(xlab="Time (Monthly)", outer=TRUE)
par(opar)
##5 point moving average
plot(T_max, type = 'o', pch = 20, ylab = "in ~(degree~C)",
     main = "Monthly average temp of Pokhara")
T_max.m5 <- filter(T_max, filter = rep(1, 5) / 5) # note NAs
lines(T_max.m5, col = "hotpink", lwd = 2)
### Spencer's 15-point filter
sp <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320
T_max_sp <- filter(T_max, filter = sp)
lines(T_max_sp, col = "steelblue", lwd = 2)
# for minimum
##5 point moving average
plot(T_min, type = 'o', pch = 20, ylab = "in ~(degree~C)",
     main = "Monthly average temp of Pokhara")
T_min.m5 <- filter(T_min, filter = rep(1, 5) / 5) # note NAs
lines(T_min.m5, col = "hotpink", lwd = 2)
### Spencer's 15-point filter
sp <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320
T_min_sp <- filter(T_min, filter = sp)
lines(T_min_sp, col = "steelblue", lwd = 2)
# residuals
#max
T_max_resid = T_max - T_max_sp
plot(T_max_resid , type = 'o', pch = 20, ylab = "Residuals",
     main = "Pokhara max Temp residuals from Spencer's 15-point moving average filter")
#min
T_min_resid = T_min - T_min_sp
plot(T_min_resid , type = 'o', pch = 20, ylab = "Residuals",
     main = "Pokhara min Temp residuals from Spencer's 15-point moving average filter")
  # check the independence of residuals
plot(c(head(T_max_resid, -1)), c(tail(T_max_resid, -1)), pch = 20,
     xlab = expression(epsilon[t]), ylab = expression(epsilon[t-1]),
     main = "Lag-1 correlation")
plot(c(head(T_max_resid, -2)), c(tail(T_max_resid, -2)), pch = 20,
     xlab = expression(epsilon[t]), ylab = expression(epsilon[t-2]),
     main = "Lag-2 correlation")
##
plot(T_max, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Monthly average temp of Pokhara")

### scatter plot smoothing default f = 2/3(smoother span propotion of points in the plot which #influence the smooth at each value. Larger values give more smoothness.)
lines(lowess(T_max), lwd = 2, col = "tan")
lines(lowess(T_max, f = 0.25), lwd = 2, col = "hotpink")
###
### filter out seasonality
f12 <- c(0.5, rep(1, 11), 0.5)
f12 <- f12 / sum(f12)
T_max_adj <- filter(T_max, filter = f12)
plot(T_max, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Monthly average temp of Pokhara")

lines(T_max_adj, col = "thistle", type = 'o', pch = 20) # deseasonalized series
T_max_ma5 <- filter(T_max_adj, filter = rep(1, 5) / 5)
lines(T_max_ma5, col = "hotpink", lwd = 2)

T_max_sp <- filter(T_max_adj, filter = sp)
lines(T_max_sp, col = "steelblue", lwd = 2)

### using built-in black-box tools
plot(stl(T_max, s.window = "periodic"))

### residuals
T_max_res <- T_max - T_max_ma5
plot(T_max_res, main = "T_max - residuals")

##Fit a structural model for a time series by maximum likelihood.
plot(T_max, xlab = "Year and Month" , ylab = expression(T[~t]~~~~(degree~C)), type = "o", pch = 16)
## Using base R functions
#The simplest model is the local level model specified by type = "level"
fitT_max <- StructTS(T_max, "level")
fitT_max
lines(fitted(fitT_max), lty = "dashed", lwd = 2)
lines(tsSmooth(fitT_max), lty = "dotted", lwd = 2)
tsdiag(fitT_max)
#The local linear trend model, type = "trend", has the same measurement equation, but with a #time-varying slope in the dynamics for μt
plot(T_max, xlab = "Year and Month" , ylab = expression(T[~t]~~~~(degree~C)), type = "o", pch = 16)
fit <- StructTS(T_max, type = "trend")
fit
plot(cbind(fitted(fit), resids = resid(fit)), main = " Pokhara max Temp in degree celcius")
tsdiag(fit)
##The basic structural model, type = "BSM", is a local trend model with an additional seasonal #component.
fit <- StructTS(T_max, type = "BSM")
plot(T_max, xlab = "Year and Month" , ylab = expression(T[~t]~~~~(degree~C)), type = "o", pch = 16)
plot(cbind(fitted(fit), resids = resid(fit)), main = " Pokhara max Temp in degree celcius")
tsdiag(fit)
#DLM
modT_max <- dlm(FF = 1, GG = 1, V = 1, W = 6.687, m0 = 1000, C0 = 1e6)
FF(modT_max) 
V(modT_max)
## Kalman Filter
filterT_max <- dlmFilter(T_max, modT_max)
str(filterT_max, 1) # everything you are getting...
## a useful function to plot filtered means and probability bounds
lines.dlmFiltered <- function(filt_obj, prob = 0.90) {
    ## INPUTS
    ## filt_obj: an object of class 'dlmFiltered', returned by 'dlmFilter'
    ## prob:     the probability level of probability bands
    ##
    ## Adds filtered values and probability bands to the current plot
    
    y <- as.ts(filt_obj$y)
    F <- filt_obj$mod$FF

    ## recover standard deviations from SVD of covariance matrices
    sdFilt <- sapply(with(filt_obj, dlmSvd2var(U.C, D.C)[-1]), function(x) F %*% x %*% t(F))
    sdFilt <- ts(sqrt(sdFilt))
    tsp(sdFilt) <- tsp(y)
    ## compute expected noise-free observation
    expected <- ts(drop(as.matrix(filt_obj$m)[-1, ] %*% t(F)))
    tsp(expected) <- tsp(y)
    
    alpha <- 0.5 * (1 - prob) # tail probability

    ## compute probability limits and plot, together with filtering mean
    lns <- cbind(expected, qnorm(alpha, sd = sdFilt) %o% c(-1, 1) + as.vector(expected)) 
    lines(lns[, 1], col = "blue", lty = "longdash", lwd = 2)
    for (i in 1 : 2) lines(lns[, i + 1], col = "red", lty = "dotdash", lwd = 2)
    invisible()
}

## try it out
plot(T_max, ylab = "Temperature", type = 'o', pch = 16)
lines(filterT_max)

## Diagnostic plots
tsdiag(filterT_max)
qqnorm(residuals(filterT_max, sd = FALSE))
qqline(residuals(filterT_max, sd = FALSE))
###
stop()
### Spencer's 15-point filter
sp <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320
huron_sp <- filter(LakeHuron, filter = sp)
lines(huron_sp, col = "steelblue", lwd = 2)

## a bit fancier
opar <- par(mfrow=c(3,1), mar=c(0,3.5,0,3), oma=c(3.5,0,2,0), mgp=c(2,.6,0), cex.lab=1.1, tcl=-.3, las=1)
plot(cmort, ylab=expression(M[~t]), xaxt="no", type='n')
  grid(lty=1, col=gray(.9))
  lines(cmort, col=rgb(0,0,.9))
plot(tempr, ylab=expression(T[~t]~~~~(degree~F)), xaxt="no", yaxt='no', type='n')
  grid(lty=1, col=gray(.9))
  lines(tempr, col=rgb(.9,0,.9)) 
  axis(4) 
plot(part, ylab=expression(P[~t]~~~~(PPM)))
  grid(lty=1, col=gray(.9))
  lines(part, col=rgb(0,.7,0))
title(xlab="Time (weekly)", outer=TRUE)
par(opar)
##
plot(LakeHuron, type = 'o', pch = 20, ylab = "Level (feet)",
     main = "Annual measurements of Lake Huron level")

str(LakeHuron) # an object of class "ts"
class(LakeHuron)
start(LakeHuron)
end(LakeHuron)
frequency(LakeHuron)
tsp(LakeHuron) # start, end, and frequency at once


huron_ma5 <- filter(LakeHuron, filter = rep(1, 5) / 5) # note NAs
lines(huron_ma5, col = "hotpink", lwd = 2)
##
## residuals
### Spencer's 15-point filter
sp <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320
huron_sp <- filter(LakeHuron, filter = sp)
lines(huron_sp, col = "steelblue", lwd = 2)


huron_resid <- LakeHuron - huron_sp
plot(huron_resid, type = 'o', pch = 20, ylab = "Residuals",
     main = "Lake Huron level\nResiduals from Spencer's 15-point moving average filter")

### Residuals are not independend...
plot(c(head(huron_resid, -1)), c(tail(huron_resid, -1)), pch = 20,
     xlab = expression(epsilon[t]), ylab = expression(epsilon[t-1]),
     main = "Lag-1 correlation")
plot(c(head(huron_resid, -2)), c(tail(huron_resid, -2)), pch = 20,
     xlab = expression(epsilon[t]), ylab = expression(epsilon[t-2]),
     main = "Lag-2 correlation")
plot(globtemp, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Global mean land-ocean temperature\nDeparture from 1951-1980 average")

### 'lowess'
lines(lowess(globtemp), lwd = 2, col = "tan")
lines(lowess(globtemp, f = 0.25), lwd = 2, col = "hotpink")

### 'supsmu'
lines(supsmu(time(globtemp), globtemp), lwd = 2, col = "steelblue")

legend("topleft", legend = c("lowess", "lowess (\'f = 0.25\')", "supsmu"),
       lwd = 2, col = c("tan", "hotpink", "steelblue"), bty = 'n',
       title = "Local smoothers")

####
plot(wineind, type = 'o', pch = 20, ylab = expression(10^3 ~ "liters"),
     main = "Australian wine sales (monthly)")

### filter out seasonality
f12 <- c(0.5, rep(1, 11), 0.5)
f12 <- f12 / sum(f12)
wine_adj <- filter(wineind, filter = f12)

lines(wine_adj, col = "thistle", type = 'o', pch = 20) # deseasonalized series
wine_ma5 <- filter(wine_adj, filter = rep(1, 5) / 5)
lines(wine_ma5, col = "hotpink", lwd = 2)

wine_sp <- filter(wine_adj, filter = sp)
lines(wine_sp, col = "steelblue", lwd = 2)

### using built-in black-box tools
plot(stl(wineind, s.window = "periodic"))

### residuals
wine_res <- wineind - wine_ma5
plot(wine_res, main = "Wine sales - residuals")




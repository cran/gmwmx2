## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 6,
  fig.align = "center"
)



## -----------------------------------------------------------------------------
library(gmwmx2)
set.seed(123)

## -----------------------------------------------------------------------------
model_wn <- wn(sigma2 = 1)
series_wn <- generate(model_wn, n = 500)
plot(series_wn)

## -----------------------------------------------------------------------------
head(series_wn$series)

## -----------------------------------------------------------------------------
model_ar1 <- ar1(phi = 0.8, sigma2 = 1)
series_ar1 <- generate(model_ar1, n = 500)
plot(series_ar1)

## -----------------------------------------------------------------------------
model_pl <- pl(kappa = 0.3, sigma2 = 1)
series_pl <- generate(model_pl, n = 500)
plot(series_pl)

## -----------------------------------------------------------------------------
model_matern <- matern(sigma2 = 1, lambda = 1, alpha=1)
series_matern <- generate(model_matern, n = 500)
plot(series_matern)

## -----------------------------------------------------------------------------
model_rw <- rw(sigma2 = 0.1)
series_rw <- generate(model_rw, n = 500)
plot(series_rw)

## -----------------------------------------------------------------------------
model_fl <- flicker(sigma2 = 1)
series_fl <- generate(model_fl, n = 500)
plot(series_fl)

## -----------------------------------------------------------------------------
model_wn <- wn(sigma2 = 1)

set.seed(1234)
series_a <- generate(model_wn, n = 100)

series_b <- generate(model_wn, n = 100, seed = 1234)
series_c <- generate(model_wn, n = 100, seed = 4321)

all.equal(series_a$series, series_b$series)
all.equal(series_a$series, series_c$series)

## -----------------------------------------------------------------------------
model_wn_ar1 <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
class(model_wn_ar1)
series_wn_ar1 <- generate(model_wn_ar1, n = 500)
plot(series_wn_ar1)


## -----------------------------------------------------------------------------
names(series_wn_ar1)


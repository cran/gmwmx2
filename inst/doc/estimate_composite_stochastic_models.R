## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 6,
  fig.align = "center"
)


knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
# library(gmwmx2)
# boxplot_mean_dot <- function(x, ...) {
#   boxplot(x, ...)
#   x_mat <- as.matrix(x)
#   mean_vals <- colMeans(x_mat, na.rm = TRUE)
#   points(seq_along(mean_vals), mean_vals, pch = 16, col = "black", cex = 1.5)
# }

## -----------------------------------------------------------------------------
# n <- 10000
# sigma2_wn = 25
# phi_ar1 = 0.99
# sigma2_ar1 = 1
# mod1 <- wn(sigma2 = sigma2_wn) + ar1(phi = phi_ar1, sigma2 = sigma2_ar1)
# y1 <- generate(mod1, n = n)
# plot(y1)
# fit1 <- gmwm2(y1, model = wn() + ar1())
# fit1
# plot(fit1)

## -----------------------------------------------------------------------------
# B = 100
# mat_res = matrix(NA, nrow = B, ncol = 3)
# for(b in seq(B)){
#   y <- generate(mod1, n = n)
#   fit = gmwm2(y, model = wn() + ar1())
#   mat_res[b,] = c(fit$theta_domain$`White Noise_1`, fit$theta_domain$`AR(1)_2`)
#   # cat("Done with b =", b, "\n")
# }

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 3))
# boxplot_mean_dot(mat_res[, 1], main = expression(sigma[WN]^2), ylab = "Estimated Value")
# abline(h = sigma2_wn)
# boxplot_mean_dot(mat_res[, 2], main = expression(phi[AR1]), ylab = "Estimated Value")
# abline(h = phi_ar1)
# boxplot_mean_dot(mat_res[, 3], main = expression(sigma[AR1]^2), ylab = "Estimated Value")
# abline(h = sigma2_ar1)
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# sigma2_wn = 5
# kappa_pl = -0.9
# sigma2_pl = 1
# mod2 <- wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl)
# y2 <- generate(mod2, n = n)
# plot(y2)
# fit2 <- gmwm2(y2, model = wn() + pl())
# fit2
# plot(fit2)

## -----------------------------------------------------------------------------
# B = 100
# mat_res = matrix(NA, nrow = B, ncol = 3)
# for(b in seq(B)){
#   y <- generate(mod2, n = n)
#   fit2 = gmwm2(y, model = wn() + pl())
#   mat_res[b,] = c(fit2$theta_domain$`White Noise_1`, fit2$theta_domain$`Stationary PowerLaw_2`)
#   # cat("Done with b =", b, "\n")
# }

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 3))
# boxplot_mean_dot(mat_res[, 1], main = expression(sigma[WN]^2), ylab = "Estimated Value")
# abline(h = sigma2_wn)
# boxplot_mean_dot(mat_res[, 2], main = expression(kappa[PL]), ylab = "Estimated Value")
# abline(h = kappa_pl)
# boxplot_mean_dot(mat_res[, 3], main = expression(sigma[PL]^2), ylab = "Estimated Value")
# abline(h = sigma2_pl)
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# sigma2_wn = 5
# phi_ar1 = 0.98
# sigma2_ar1 = 1
# sigma2_rw = 0.1
# mod3 <- wn(sigma2_wn) + ar1(phi = phi_ar1, sigma2 = sigma2_ar1) + rw(sigma2_rw)
# y3 <- generate(mod3, n = n)
# plot(y3)
# fit3 <- gmwm2(y3, model = wn() + ar1() + rw())
# fit3
# plot(fit3)

## -----------------------------------------------------------------------------
# B = 100
# mat_res = matrix(NA, nrow = B, ncol = 4)
# for(b in seq(B)){
#   y <- generate(mod3, n = n)
#   fit3 = gmwm2(y, model = wn() + ar1() + rw())
#   mat_res[b,] = c(fit3$theta_domain$`White Noise_1`, fit3$theta_domain$`AR(1)_2`, fit3$theta_domain$`Random Walk_3`)
#   # cat("Done with b =", b, "\n")
# }

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 4))
# boxplot_mean_dot(mat_res[, 1], main = expression(sigma[WN]^2), ylab = "Estimated Value")
# abline(h = sigma2_wn)
# boxplot_mean_dot(mat_res[, 2], main = expression(phi[AR1]), ylab = "Estimated Value")
# abline(h = phi_ar1)
# boxplot_mean_dot(mat_res[, 3], main = expression(sigma[AR1]^2), ylab = "Estimated Value")
# abline(h = sigma2_ar1)
# boxplot_mean_dot(mat_res[, 4], main = expression(sigma[RW]^2), ylab = "Estimated Value")
# abline(h = sigma2_rw)
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# sigma2_wn = 20
# sigma2_rw = 0.1
# mod4 <- wn(sigma2_wn) + rw(sigma2_rw)
# y4 <- generate(mod4, n = n)
# plot(y4)
# fit4 <- gmwm2(y4, model = wn() + rw())
# fit4
# plot(fit4)

## -----------------------------------------------------------------------------
# B = 100
# mat_res = matrix(NA, nrow = B, ncol = 2)
# for(b in seq(B)){
#   y <- generate(mod4, n = n)
#   fit4 = gmwm2(y, model = wn()  + rw())
#   mat_res[b,] = c(fit4$theta_domain$`White Noise_1`, fit4$theta_domain$`Random Walk_2`)
#   # cat("Done with b =", b, "\n")
# }

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 2))
# boxplot_mean_dot(mat_res[, 1], main = expression(sigma[WN]^2), ylab = "Estimated Value")
# abline(h = sigma2_wn)
# boxplot_mean_dot(mat_res[, 2], main = expression(sigma[RW]^2), ylab = "Estimated Value")
# abline(h = sigma2_rw)
# par(mfrow = c(1, 1))


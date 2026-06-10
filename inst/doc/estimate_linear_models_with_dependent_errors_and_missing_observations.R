## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 6,
  fig.align = "center"
)

# so as not to execute this vignette during package development
knitr::opts_chunk$set(eval = FALSE)

## ----warning=FALSE, message=FALSE---------------------------------------------
# library(gmwmx2)
# library(wv)
# library(dplyr)
# 
# boxplot_mean_dot <- function(x, ...) {
#   boxplot(x, ...)
#   x_mat <- as.matrix(x)
#   mean_vals <- colMeans(x_mat, na.rm = TRUE)
#   points(seq_along(mean_vals), mean_vals, pch = 16, col = "black")
# }

## -----------------------------------------------------------------------------
# n = 15*365
# X = matrix(NA, nrow = n, ncol = 4)
# # intercept
# X[, 1] = 1
# # trend
# X[, 2] = 1:n
# # annual sinusoid
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(0.5 , 0.00004, 0.005,  0.0008)
# 
# # visualize the deterministic signal
# plot(x = X[, 2], y = X %*% beta, type = "l",
#      main = "Signal", xlab = "Time", ylab = "Signal")

## -----------------------------------------------------------------------------
# phi_ar1 = 0.975
# sigma2_ar1 =  5e-05
# sigma2_wn =  7e-04
# eps = generate(ar1(phi = phi_ar1, sigma2 = sigma2_ar1) + wn(sigma2_wn), n = n, seed = 123)$series
# plot(wv::wvar(eps))
# y = X %*% beta + eps
# 
# # generate missingness
# mod_missing = markov_two_states(p1 = 0.05, p2 = .95)
# Z = generate(mod_missing, n = n, seed = 123)
# plot(Z)
# Z_process = Z$series
# y_observed = ifelse(Z_process == 1, y, NA)
# plot(X[, 2], y_observed, type = "l", main = "Simulated y (WN + AR1)")
# 
# 
# 
# id_not_observed = which(is.na(y_observed))
# for(i in id_not_observed){
#   abline(v = X[i, 2], col = "grey70", lty =1)
# }

## -----------------------------------------------------------------------------
# fit = gmwmx2(X = X, y = y_observed, model = wn() + ar1() )
# print(fit)

## -----------------------------------------------------------------------------
# B = 100
# mat_res = matrix(NA, nrow=B, ncol=19)
# for(b in seq(B)){
# 
#   eps = generate(ar1(phi=phi_ar1, sigma2=sigma2_ar1) + wn(sigma2_wn), n=n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   Z = generate(mod_missing, n = n, seed = 123 + b)$series
#   y_observed = ifelse(Z == 1, y, NA)
#   fit = gmwmx2(X = X, y = y_observed, model = wn() + ar1() )
#   # misspecified model assuming white noise as the stochastic model and only on observed data
#   fit2 = lm(y_observed~X[,2] + X[,3] + X[,4])
# 
#   mat_res[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                    summary(fit2)$coefficients[,1],
#                    summary(fit2)$coefficients[,2],
#                    fit$theta_domain$`AR(1)_2`,
#                    fit$theta_domain$`White Noise_1`)
#   # cat("Iteration ", b, " completed.\n")
# }

## -----------------------------------------------------------------------------
# # compute empirical coverage
# mat_res_df = as.data.frame(mat_res)
# colnames(mat_res_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                          "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                          "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                          "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                          "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                          "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                          "phi_ar1","sigma_2_ar1" ,"sigma_2_wn")
# 
# 
# 
# par(mfrow = c(1, 4))
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#                  names = c("beta0"),
#                  main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#                  names = c("beta1"),
#                  main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#                  names = c("beta2"),
#                  main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#                  names = c("beta3"),
#                  main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
# 
# 
# par(mfrow = c(1, 3))
# 
# boxplot_mean_dot(mat_res_df$phi_ar1,las=1,
#                  names = c("phi_ar1"),
#                  main = expression(phi["AR1"]), ylab = "Estimate")
# abline(h = phi_ar1, col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df$sigma_2_ar1,las=1,
#                  names = c("phi_ar1"),
#                  main = expression(sigma["AR1"]^2), ylab = "Estimate")
# abline(h = sigma2_ar1, col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df$sigma_2_wn,las=1,
#                  names = c("phi_ar1"),
#                  main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn, col = "black", lwd = 2)
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# # Compute empirical coverage of confidence intervals for beta
# zval = qnorm(0.975)
# mat_res_df$upper_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat + zval * mat_res_df$gmwmx_std_beta1_hat
# mat_res_df$lower_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat - zval * mat_res_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[2], B), mat_res_df$lower_ci_gmwmx_beta1, mat_res_df$upper_ci_gmwmx_beta1) %>% mean()
# 
# # do the same for lm beta
# mat_res_df$upper_ci_lm_beta1 = mat_res_df$lm_beta1_hat + zval * mat_res_df$lm_std_beta1_hat
# mat_res_df$lower_ci_lm_beta1 = mat_res_df$lm_beta1_hat - zval * mat_res_df$lm_std_beta1_hat
# dplyr::between(rep(beta[2], B), mat_res_df$lower_ci_lm_beta1, mat_res_df$upper_ci_lm_beta1) %>% mean()

## -----------------------------------------------------------------------------
# kappa_pl <- -0.75
# sigma2_wn <- 2e-07
# sigma2_pl <- 2.25e-06
# eps_pl = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl),
#                   n = n, seed = 1234)
# plot(eps_pl$series, type = "l")
# plot(wv::wvar(eps_pl$series))
# y_pl = X %*% beta + eps_pl$series
# Z = generate(mod_missing, n = n, seed = 1234)
# plot(Z)
# Z_process = Z$series
# y_observed = ifelse(Z_process == 1, y_pl, NA)
# 
# plot(X[, 2], y_observed, type = "l", main = "Simulated y (WN + PL)")
# id_not_observed = which(is.na(y_observed))
# for(i in id_not_observed){
#   abline(v = X[i, 2], col = "grey70", lty =1)
# }

## -----------------------------------------------------------------------------
# fit_pl = gmwmx2(X = X, y = y_observed, model = wn() + pl())
# fit_pl

## -----------------------------------------------------------------------------
# B_pl = 100
# mat_res_pl = matrix(NA, nrow = B_pl, ncol = 19)
# for(b in seq(B_pl)){
#   eps = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl),
#                  n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
# 
#   Z = generate(mod_missing, n = n, seed = (123 + b))
#   Z_process = Z$series
#   y_observed = ifelse(Z_process == 1, y, NA)
# 
#   fit = gmwmx2(X = X, y = y_observed, model = wn() + pl())
#   fit2 = lm(y_observed~X[,2] + X[,3] + X[,4])
# 
#   mat_res_pl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Stationary PowerLaw_2`,
#                       fit$theta_domain$`White Noise_1`)
#   # cat("Iteration ", b, " \n")
# }

## -----------------------------------------------------------------------------
# mat_res_pl_df = as.data.frame(mat_res_pl)
# colnames(mat_res_pl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat",
#                             "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "kappa_pl","sigma2_pl" ,"sigma2_wn")
# 
# 
# # Plot empirical distributions (beta and stochastic parameters)
# 
# 
# par(mfrow = c(1, 4))
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#                  names = c("beta0"),
#                  main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#                  names = c("beta1"),
#                  main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#                  names = c("beta2"),
#                  main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#                  names = c("beta3"),
#                  main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
# 
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 3))
# 
# boxplot_mean_dot(mat_res_pl_df$kappa_pl,las=1,
#                  main = expression(kappa["PL"]), ylab = "Estimate")
# abline(h = kappa_pl, col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_pl_df$sigma2_pl,las=1,
#                  main = expression(sigma["PL"]^2), ylab = "Estimate")
# abline(h = sigma2_pl, col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_pl_df$sigma2_wn,las=1,
#                  names = c("phi_ar1"),
#                  main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn, col = "black", lwd = 2)
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# zval = qnorm(0.975)
# mat_res_pl_df$upper_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat + zval * mat_res_pl_df$gmwmx_std_beta1_hat
# mat_res_pl_df$lower_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat - zval * mat_res_pl_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[2], B_pl), mat_res_pl_df$lower_ci_gmwmx_beta1, mat_res_pl_df$upper_ci_gmwmx_beta1) %>% mean()
# 
# # do the same for lm beta
# mat_res_pl_df$upper_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat + zval * mat_res_pl_df$lm_std_beta1_hat
# mat_res_pl_df$lower_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat - zval * mat_res_pl_df$lm_std_beta1_hat
# 
# dplyr::between(rep(beta[2], B_pl), mat_res_pl_df$lower_ci_lm_beta1, mat_res_pl_df$upper_ci_lm_beta1) %>% mean()

## -----------------------------------------------------------------------------
# sigma2_wn_fl <- 8e-07
# sigma2_fl <- 1e-06
# eps_fl = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl),
#                   n = n, seed = 123)
# plot(eps_fl$series, type = "l")
# plot(wv::wvar(eps_fl$series))
# Z = generate(mod_missing, n = n, seed = 1234)
# plot(Z)
# 
# y_fl = X %*% beta + eps_fl$series
# 
# Z = generate(mod_missing, n = n, seed = 1234)
# plot(Z)
# Z_process = Z$series
# 
# 
# y_observed = ifelse(Z_process == 1, y_fl, NA)

## -----------------------------------------------------------------------------
# plot(X[, 2], y_observed, type = "l", main = "Simulated y (WN + FL)")
# id_not_observed = which(is.na(y_observed))
# for(i in id_not_observed){
#   abline(v = X[i, 2], col = "grey70", lty =1)
# }

## -----------------------------------------------------------------------------
# fit_fl = gmwmx2(X = X, y = y_observed, model = wn() + flicker())
# fit_fl

## -----------------------------------------------------------------------------
# B_fl = 100
# mat_res_fl = matrix(NA, nrow = B_fl, ncol = 18)
# for(b in seq(B_fl)){
#   eps = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl),
#                  n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
# 
# 
#   Z = generate(mod_missing, n = n, seed = (123 + b))
#   Z_process = Z$series
#   y_observed = ifelse(Z_process == 1, y, NA)
# 
# 
# 
#   fit = gmwmx2(X = X, y = y_observed, model = wn() + flicker())
#   fit2 = lm(y_observed~X[,2] + X[,3] + X[,4])
# 
#   mat_res_fl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Flicker_2`,
#                       fit$theta_domain$`White Noise_1`)
#   # cat("Iteration ", b, " \n")
# }

## -----------------------------------------------------------------------------
# mat_res_fl_df = as.data.frame(mat_res_fl)
# colnames(mat_res_fl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat",
#                             "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "sigma2_fl" ,"sigma2_wn")

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 4))
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#                  names = c("beta0"),
#                  main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#                  names = c("beta1"),
#                  main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#                  names = c("beta2"),
#                  main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#                  names = c("beta3"),
#                  main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
# 
# par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
# par(mfrow = c(1, 2))
# 
# boxplot_mean_dot(mat_res_fl_df$sigma2_fl,las=1,
#                  main = expression(sigma["FL"]^2), ylab = "Estimate")
# abline(h = sigma2_fl, col = "black", lwd = 2)
# 
# boxplot_mean_dot(mat_res_fl_df$sigma2_wn,las=1,
#                  names = c("phi_ar1"),
#                  main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn_fl, col = "black", lwd = 2)
# par(mfrow = c(1, 1))
# 

## -----------------------------------------------------------------------------
# # Compute empirical coverage of confidence intervals for beta
# 
# zval = qnorm(0.975)
# mat_res_fl_df$upper_ci_gmwmx_beta1 = mat_res_fl_df$gmwmx_beta1_hat + zval * mat_res_fl_df$gmwmx_std_beta1_hat
# mat_res_fl_df$lower_ci_gmwmx_beta1 = mat_res_fl_df$gmwmx_beta1_hat - zval * mat_res_fl_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[2], B_fl), mat_res_fl_df$lower_ci_gmwmx_beta1, mat_res_fl_df$upper_ci_gmwmx_beta1) %>% mean()
# 
# # do the same for lm beta
# mat_res_fl_df$upper_ci_lm_beta1 = mat_res_fl_df$lm_beta1_hat + zval * mat_res_fl_df$lm_std_beta1_hat
# mat_res_fl_df$lower_ci_lm_beta1 = mat_res_fl_df$lm_beta1_hat - zval * mat_res_fl_df$lm_std_beta1_hat
# dplyr::between(rep(beta[2], B_fl), mat_res_fl_df$lower_ci_lm_beta1, mat_res_fl_df$upper_ci_lm_beta1) %>% mean()
# 


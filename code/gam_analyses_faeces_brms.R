library(data.table)
library(mgcv)
library(gratia)
library(brms)
library(rstan)

# options
rstan_options(auto_write = TRUE)
options(mc.cores = (parallel::detectCores() / 2))

if(basename(getwd())!="typhoid_data") setwd("~/Documents/sf_drugs_projects/typhoid_data")
while (!is.null(dev.list())) dev.off()
rm(list=ls())

## Source functions
source("~/Documents/sf_drugs_projects/typhoid_data/code/gam_analyses_fn.R")

## Load data
use.imputed <- TRUE
if (!use.imputed) {
  load("./data/formatted_datatables.RData")
} else {
  load("./data/formatted_datatables_imputedwt.RData")
}
data.Ab <- copy(data.O) # insert the used data frame here
## Cro, Ce, A have insufficient data for GAM

## Additional variables
data.Ab[,extraDur:=(DurAb1 - DurRand)]
data.Ab[,doseperMIC:=dose/OfMIC]
data.Ab[,Age.scaled:=scale(Age)]

## remove unused factors of sty1
data.Ab$sty1 <- factor(data.Ab$sty1)

## 1. outcome: FCT38
## Check difference between missing and included data
data.Ab[,.(
  mean(totaldose_alloc,na.rm=TRUE),
  sd(totaldose_alloc,na.rm=TRUE) / sqrt(sum(!is.na(totaldose_alloc)))
), by = is.na(F2P)]
data.Ab[,.(mean(OfMIC,na.rm=TRUE),sd(OfMIC,na.rm=TRUE) / sqrt(sum(!is.na(OfMIC)))),by = is.na(F2P)]
data.Ab[,.(mean(wt,na.rm=TRUE),sd(wt,na.rm=TRUE) / sqrt(sum(!is.na(wt)))),by = is.na(F2P)]

## We first want to account for the fact that we know a bit more than this
## - some are pos, but we don't know how many of their samples were pos
## First create a test model, to get the stan code and data
data.Ab.brms <- data.Ab[!is.na(F2P) & !is.na(OfMIC) & !is.na(totaldose_alloc)]
data.Ab.brms[F2P==0, F2R:=0]
data.Ab.brms[F2P==1 & is.na(F2R), F2R:=1] ## Make 1 so stan will generate code; edit later 

# create a brms model, but only for the data and object - does not need to be fit
# (do this rather than make_standata, as need the object for prediction on newdata
gam.Ab.brms <- brm(
  F2R | trials(F2N) ~ t2(totaldose_alloc, OfMIC) + s(wt) + s(sty1,bs="re"),
  data = data.Ab.brms, family = binomial,
  empty = TRUE
)

# create stan data (without compiling model
data.Ab.stan <- standata(gam.Ab.brms)
# code
code.Ab.stan <- stancode(gam.Ab.brms)

# now edit the likelihood in this file to account for when there is at least one positive
# using binomial_lccdf() (log complementary cumulative distribution function)
code.Ab.stan <- gsub(
  "\\Qtarget += binomial_logit_lpmf(Y | trials, mu);\\E",
  "for (ii in 1:N) {
      if (Y[ii] == -1) {
        target += binomial_lccdf(0 | trials[ii], inv_logit(mu[ii]));
      } else {
        target += binomial_logit_lpmf(Y[ii] | trials[ii], mu[ii]);
      }
    }",
  code.Ab.stan,
  perl = TRUE
)
code.Ab.stan <- gsub(
  "\\Q\n  if (!prior_only) {\n    // initialize linear predictor term\n    vector[N] mu = rep_vector(0.0, N);\n    mu += Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1;\\E",
  "\n  if (!prior_only) {",
  code.Ab.stan,
  perl = TRUE
)
code.Ab.stan <- gsub(
  "\\Q}\nmodel {\\E",
  "  // initialize linear predictor term\n  vector[N] mu = rep_vector(0.0, N);\n  mu += Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1;\n}\nmodel {",
  code.Ab.stan,
  perl = TRUE
)
fc <- file("code/faeces_gam_stan.stan")
writeLines(code.Ab.stan, fc)
close(fc)

## Now fit the full model, first recreating the data and making it -1 where missing num pos
data.Ab.stanprep <- data.Ab[!is.na(F2P) & !is.na(OfMIC) & !is.na(totaldose_alloc)]
data.Ab.stanprep[F2P==0, F2R:=0]
data.Ab.stanprep[F2P==1 & is.na(F2R), F2R:=-1]
data.Ab.stan$Y <- data.Ab.stanprep$F2R

niter <- 4000
nwarmup <- 2000
nchains <- 4
gam.Ab.faeces <- stan(
  "code/faeces_gam_stan.stan",
  iter = niter,
  warmup = nwarmup,
  data = data.Ab.stan,
  chains = nchains,
  cores = nchains,
  seed = 161,
  refresh = 100,
  control = list(adapt_delta = 0.99)
)

# hmc diagnostics
check_hmc_diagnostics(gam.Ab.faeces)

# traceplots
traceplot(gam.Ab.faeces, pars = "Intercept")
traceplot(gam.Ab.faeces, pars = "bs")
traceplot(gam.Ab.faeces, pars = "zs_1_1")
traceplot(gam.Ab.faeces, pars = "zs_1_2")
traceplot(gam.Ab.faeces, pars = "zs_1_3")
traceplot(gam.Ab.faeces, pars = "zs_1_1")
traceplot(gam.Ab.faeces, pars = "zs_2_1")
traceplot(gam.Ab.faeces, pars = "zs_3_1")
traceplot(gam.Ab.faeces, pars = "sds_1")
traceplot(gam.Ab.faeces, pars = "sds_2")
traceplot(gam.Ab.faeces, pars = "sds_3")
traceplot(gam.Ab.faeces, pars = "s_1_1")
traceplot(gam.Ab.faeces, pars = "s_1_2")
traceplot(gam.Ab.faeces, pars = "s_1_3")
traceplot(gam.Ab.faeces, pars = "s_2_1")
traceplot(gam.Ab.faeces, pars = "s_3_1")

gam.Ab.faeces_draws <- as_draws_df(gam.Ab.faeces)

## Create posterior g-formula average treatment effect
dose.vec <- seq(
  floor(min(data.Ab.stanprep$totaldose_alloc)),
  ceiling(max(data.Ab.stanprep$totaldose_alloc)),
  1
)
# set up data structures
quantiles_desired <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
mu_ate <- data.frame(matrix(
  NA, nrow = length(dose.vec), ncol = 1 + length(quantiles_desired)
))
colnames(mu_ate) <- c("mean", quantiles_desired)
p_ate <- post_ate <- mu_ate
# create baseline pred at mean dose (53.7)
Edose <- mean(data.Ab.stanprep$totaldose_alloc, na.rm=TRUE)
data.Ab.temp <- copy(data.Ab.brms)
data.Ab.temp$totaldose_alloc <- Edose
data.Ab.stan.temp <- standata(gam.Ab.brms, newdata = data.Ab.temp)
post.temp <- create_posterior_pred_faeces(data.Ab.stan.temp, gam.Ab.faeces_draws)
mu_base <- post.temp$mu
p_base <- post.temp$p
post_base <- post.temp$post
# now try for the range of counterfactual doses
data.Ab.temp <- copy(data.Ab.brms)
for (ii in seq_along(dose.vec)) {
  print(ii)
  data.Ab.temp$totaldose_alloc <- dose.vec[ii]
  data.Ab.stan.temp <- standata(gam.Ab.brms, newdata = data.Ab.temp)
  post.temp <- create_posterior_pred_faeces(data.Ab.stan.temp, gam.Ab.faeces_draws)
  mu_ate.temp <- rowMeans(post.temp$mu - mu_base)
  mu_ate[ii, 1] <- mean(mu_ate.temp)
  mu_ate[ii, -1] <- quantile(mu_ate.temp, probs = quantiles_desired)
  p_ate.temp <- rowMeans(post.temp$p - p_base)
  p_ate[ii, 1] <- mean(p_ate.temp)
  p_ate[ii, -1] <- quantile(p_ate.temp, probs = quantiles_desired)
  post_ate.temp <- rowMeans(post.temp$post - post_base)
  post_ate[ii, 1] <- mean(post_ate.temp)
  post_ate[ii, -1] <- quantile(post_ate.temp, probs = quantiles_desired)
}

pdf("./figs/gformula_effect_faeces_brms.pdf", width = 12 / 2.54, height = 12 / 2.54, pointsize = 13)
plot(
  dose.vec,
  p_ate[, "mean"],
  ylim = range(p_ate[, c("0.05", "0.95")]),
  lty = 1,
  lwd = 2.5,
  type = "l",
  bty = "n",
  xaxs = "i",
  yaxs = "i",
  xlab = "total dose (mg/kg)",
  ylab = "Expected change in probability of positive stool sample",
  )
lines(dose.vec, p_ate[, "0.05"], lty = 3)
lines(dose.vec, p_ate[, "0.95"], lty = 3)
abline(h = 0, lty = 3)
polygon(
  c(dose.vec, rev(dose.vec)),
  c(p_ate[, "0.05"], rev(p_ate[, "0.95"])),
  col = adjustcolor("grey",0.3), border=FALSE
)
abline(h = 0, lty = 2)
dev.off()

# now do this for MIC too!!
## Create posterior g-formula average treatment effect
MIC.vec <- 2^c(-3, -2, 0)
# set up data structures
mu_ate_arr <- array(
  NA,
  dim = c(length(MIC.vec), length(dose.vec), 1 + length(quantiles_desired)),
  dimnames = list(MIC.vec, dose.vec, c("mean", quantiles_desired))
)
p_ate_arr <- post_ate_arr <- mu_ate_arr
# now try for the range of counterfactual doses
data.Ab.temp <- copy(data.Ab.brms)
for (jj in seq_along(MIC.vec)) {
  data.Ab.temp$OfMIC <- MIC.vec[jj]
  # create baseline pred at mean dose (53.7) and MIC of interest
  data.Ab.temp$totaldose_alloc <- Edose
  data.Ab.stan.temp <- standata(gam.Ab.brms, newdata = data.Ab.temp)
  post.temp <- create_posterior_pred_faeces(data.Ab.stan.temp, gam.Ab.faeces_draws)
  mu_base <- post.temp$mu
  p_base <- post.temp$p
  post_base <- post.temp$post
  for (ii in seq_along(dose.vec)) {
    print(ii)
    data.Ab.temp$totaldose_alloc <- dose.vec[ii]
    data.Ab.stan.temp <- standata(gam.Ab.brms, newdata = data.Ab.temp)
    post.temp <- create_posterior_pred_faeces(data.Ab.stan.temp, gam.Ab.faeces_draws)
    mu_ate.temp <- rowMeans(post.temp$mu - mu_base)
    mu_ate_arr[jj, ii, 1] <- mean(mu_ate.temp)
    mu_ate_arr[jj, ii, -1] <- quantile(mu_ate.temp, probs = quantiles_desired)
    p_ate.temp <- rowMeans(post.temp$p - p_base)
    p_ate_arr[jj, ii, 1] <- mean(p_ate.temp)
    p_ate_arr[jj, ii, -1] <- quantile(p_ate.temp, probs = quantiles_desired)
    post_ate.temp <- rowMeans(post.temp$post - post_base)
    post_ate_arr[jj, ii, 1] <- mean(post_ate.temp)
    post_ate_arr[jj, ii, -1] <- quantile(post_ate.temp, probs = quantiles_desired)
  }
}

pdf(
  "./figs/gformula_effect_faeces_bymic_brms.pdf",
  width = 17 / 2.54, height = 17 / 3 / 2.54, pointsize = 13
)
par(mfrow = c(1, 3), mar = c(4.1, 5.1, 2.1, 1.1))
plot(
  dose.vec,
  p_ate_arr[1, , "mean"],
  ylim = range(p_ate_arr[, , c("0.05", "0.95")]),
  lty = 1,
  lwd = 2.5,
  type = "l",
  bty = "n",
  xaxs = "i",
  yaxs = "i",
  xlab = "total dose (mg/kg)",
  ylab = "Expected change in probability\nstool sample is positive",
  main = "MIC: 0.125 mg/L"
)
lines(dose.vec, p_ate_arr[1, , "0.05"], lty = 3)
lines(dose.vec, p_ate_arr[1, , "0.95"], lty = 3)
polygon(
  c(dose.vec, rev(dose.vec)),
  c(p_ate_arr[1, , "0.05"], rev(p_ate_arr[1, , "0.95"])),
  col = adjustcolor("grey",0.3), border=FALSE
)
abline(h = 0, lty = 2)
abline(v = mean(data.Ab.brms$totaldose_alloc, na.rm = TRUE), lty = 2)
plot(
  dose.vec,
  p_ate_arr[2, , "mean"],
  ylim = range(p_ate_arr[, , c("0.05", "0.95")]),
  lty = 1,
  lwd = 2.5,
  type = "l",
  bty = "n",
  xaxs = "i",
  yaxs = "i",
  xlab = "total dose (mg/kg)",
  ylab = "",
  main = "MIC: 0.25 mg/L"
)
lines(dose.vec, p_ate_arr[2, , "0.05"], lty = 3)
lines(dose.vec, p_ate_arr[2, , "0.95"], lty = 3)
polygon(
  c(dose.vec, rev(dose.vec)),
  c(p_ate_arr[2, , "0.05"], rev(p_ate_arr[2, , "0.95"])),
  col = adjustcolor("grey",0.3), border=FALSE
)
abline(h = 0, lty = 2)
abline(v = mean(data.Ab.brms$totaldose_alloc, na.rm = TRUE), lty = 2)
plot(
  dose.vec,
  p_ate_arr[3, , "mean"],
  ylim = range(p_ate_arr[, , c("0.05", "0.95")]),
  lty = 1,
  lwd = 2.5,
  type = "l",
  bty = "n",
  xaxs = "i",
  yaxs = "i",
  xlab = "total dose (mg/kg)",
  ylab = "",
  main = "MIC: 1 mg/L"
)
lines(dose.vec, p_ate_arr[3, , "0.05"], lty = 3)
lines(dose.vec, p_ate_arr[3, , "0.95"], lty = 3)
polygon(
  c(dose.vec, rev(dose.vec)),
  c(p_ate_arr[3, , "0.05"], rev(p_ate_arr[3, , "0.95"])),
  col = adjustcolor("grey",0.3), border=FALSE
)
abline(h = 0, lty = 2)
abline(v = mean(data.Ab.brms$totaldose_alloc, na.rm = TRUE), lty = 2)
dev.off()

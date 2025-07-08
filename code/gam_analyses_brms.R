library(data.table)
library(doParallel)
library(abind)
library(brms)

# options
options(mc.cores = (parallel::detectCores() / 2))

# Clear workspace
rm(list=ls())

## Source functions
source("./code/gam_analyses_fn.R")

## Load data
load("./data/synthetic_datatables.RData")
data.Ab$sty1 <- factor(data.Ab$sty1)
  
# Remove unused columns to save space later
data.Ab <- data.Ab[, .(fct38, failure, totaldose_alloc, OfMIC, wt, sty1)]

## Pre-treatment values (Stuart & Pullen)
fct.pt <- 31.5 * 24
relapse.pt <- 45 / 360

## data.Ab$fct38 <- data.Ab$fct38 + data.Ab$Dur * 24

## 1. outcome: FCT38
## main: weight, study RE, duration, MIC, Age
get_prior(
  fct38 ~ t2(totaldose_alloc, OfMIC, k=c(NULL, NULL)) + s(wt) + s(sty1, bs="re"),
  data = data.Ab,
  family = "Gamma"
)
gam.Ab.fct.prior <- prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") +
  prior("gamma(2, 2)", class = "shape")

# get complete cases
data.Ab.fct <- data.Ab[!is.na(fct38)][!is.na(totaldose_alloc)][!is.na(OfMIC)][!is.na(wt)][!is.na(sty1)]

# fit model
gam.Ab.fct <- brm(
  fct38 ~ t2(totaldose_alloc, OfMIC, k=c(NULL, NULL)) + s(wt) + s(sty1, bs="re"),
  data = data.Ab.fct,
  family = "Gamma",
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 161,
  prior = gam.Ab.fct.prior,
  # sample_prior = "only"
  control = list(adapt_delta = 0.99)  
)


## check posterior
pp_check(gam.Ab.fct)

## Get values to report in paper
data.temp <- copy(data.Ab.fct)
data.temp[["OfMIC"]] <- 1
data.temp[["totaldose_alloc"]] <- mean(data.Ab.fct[["totaldose_alloc"]], na.rm=TRUE)
mu_base <- posterior_epred(gam.Ab.fct, newdata = data.temp)
data.temp[["totaldose_alloc"]] <- mean(data.Ab.fct[["totaldose_alloc"]], na.rm=TRUE) * 2
mu_temp <- posterior_epred(gam.Ab.fct, newdata = data.temp)
mu_ate.temp <- rowMeans(mu_temp - mu_base)
print(mean(mu_ate.temp))
print(quantile(mu_ate.temp, probs = c(0.05, 0.95)))

## plot g-computation
gcomp_list <- calculate_gcomp_brms(gam.Ab.fct, data.Ab.fct)
save(gcomp_list, file = "./data/gcomp_list_fct_brms.RData")
# load("./data/gcomp_list_fct_brms.RData")
plot_gcomp_brms(
  "./figs/gformula_effect_fct_brms.pdf", gcomp_list,
  Edose = mean(data.Ab.fct[["totaldose_alloc"]], na.rm=TRUE)
)
gcomp_list <- calculate_gcomp_brms(gam.Ab.fct, data.Ab.fct, byperc = TRUE)
plot_gcomp_brms(
  "./figs/gformula_effect_fct_byperc_brms.pdf", gcomp_list,
  Edose = 100, xlab = "% API"
)


## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(gam.Ab.fct, data.Ab.fct, bymic = TRUE, refresh = 10)
save(gcomp_list, file = "./data/gcomp_list_fct_bymic_brms.RData")
# load(file = "./data/gcomp_list_fct_bymic_brms.RData")
plot_gcomp_brms(
  "./figs/gformula_effect_fct_bymic_brms.pdf", gcomp_list,
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)
print(approx(
  gcomp_list$dose.vec, gcomp_list$mu_ate[3, , 1],
  mean(data.Ab.fct[["totaldose_alloc"]], na.rm = TRUE) * 0.5
)$y)
print(approx(
  gcomp_list$dose.vec, gcomp_list$mu_ate[3, , 2],
  mean(data.Ab.fct[["totaldose_alloc"]], na.rm = TRUE) * 0.5
)$y)
print(approx(
  gcomp_list$dose.vec, gcomp_list$mu_ate[3, , 3],
  mean(data.Ab.fct[["totaldose_alloc"]], na.rm = TRUE) * 0.5
)$y)

plot_gcomp_brms_escmid(
  "./figs/gformula_effect_fct_bymic_brms_ESCMID.pdf", gcomp_list,
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 37.22 / 2.54, height = 16 / 2.54, pointsize = 32
)


## 2. outcome: failure
## main: weight, study RE, duration, MIC, Age
get_prior(
  failure ~ t2(totaldose_alloc, OfMIC, k=c(NULL, NULL)) + s(wt) + s(sty1, bs="re"),
  data = data.Ab,
  family = "binomial"
)
gam.Ab.failure.prior <- prior("normal(0, 0.5)", class = "b") 
# get complete cases
data.Ab.failure <- data.Ab[!is.na(failure)][!is.na(totaldose_alloc)][
  !is.na(OfMIC)][!is.na(wt)][!is.na(sty1)]
data.Ab.failure$failure <- as.numeric(data.Ab.failure$failure)

# fit model
gam.Ab.failure <- brm(
  failure ~ t2(totaldose_alloc, OfMIC, k=c(NULL, NULL)) + s(wt) + s(sty1, bs="re"),
  data = data.Ab.failure,
  family = "bernoulli",
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 161,
  prior = gam.Ab.failure.prior,
  # sample_prior = "only"
  control = list(adapt_delta = 0.99)  
)

# save models
save(gam.Ab.fct, gam.Ab.failure, file = "./data/baseline_brms_models.RData")
# load(file = "./data/baseline_brms_models.RData")

## check posterior
pp_check(gam.Ab.failure)

## plot g-computation
gcomp_list <- calculate_gcomp_brms(gam.Ab.failure, data.Ab.failure)
save(gcomp_list, file = "./data/gcomp_list_failure_brms.RData")
plot_gcomp_brms(
  "./figs/gformula_effect_failure_brms.pdf", gcomp_list,
  ylab = "Expected change in failure probability"
)

## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(gam.Ab.failure, data.Ab.failure, bymic = TRUE, refresh = 30)
save(gcomp_list, file = "./data/gcomp_list_failure_bymic_brms.RData")
load(file = "./data/gcomp_list_failure_bymic_brms.RData")
plot_gcomp_brms(
  "./figs/gformula_effect_failure_bymic_brms.pdf", gcomp_list,
  ylab = "Expected change in\nfailure probability",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## 3. Compare to distributions
## Tabernero data
load("./data/tabernero_data.RData")
Ab.mean <- tabernero.individual[antibiotic=="OFX"]$mean / 100
Ab.sd <- tabernero.individual[antibiotic=="OFX"]$sd / 100
Ab.sdlog <- sqrt(log(1 + (Ab.sd / Ab.mean)^2))
Ab.meanlog <- log(Ab.mean) - Ab.sdlog^2 / 2
## Alternative distribution for Ofloxacin, with low mean
Ab.mean2 <- 93.6 / 100## Based on Cipro
Ab.sd2 <- 13.3 / 100
Ab.sdlog2 <- sqrt(log(1 + (Ab.sd2 / Ab.mean2)^2))
Ab.meanlog2 <- log(Ab.mean2) - Ab.sdlog2^2 / 2

# doses and MICs
totaldose_alloc.mean <- mean(data.Ab$totaldose_alloc, na.rm = TRUE)
papi.range <- seq(0,200) / 100
MIC.range <- 2^seq(-6, 0, 0.1)
targetdose.vec <- c(totaldose_alloc.mean / 2, totaldose_alloc.mean, totaldose_alloc.mean * 2)

## Plot tabernero and toy distribution
pdf("figs/tabernero_plus_toy_dists.pdf", width = 5, height = 2)
par(mar = c(4.1, 4.1, 1.1, 1.1))
plot(
  papi.range, dlnorm(papi.range, Ab.meanlog, Ab.sdlog),
  type = "l", lwd = 2,
  xlab = "% API", ylab = "density",
  xaxs = "i", yaxs = "i", las = 1, xlim = c(0.5, 1.5),
  bty = "n"
)
lines(papi.range, dlnorm(papi.range, Ab.meanlog2, Ab.sdlog2), lwd = 2, lty = 2, col = "red")
abline(v = c(0.90, 1.10), lty = 3)
dev.off()


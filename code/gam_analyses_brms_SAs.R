library(data.table)
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
data.Ab[, doseperMIC := dose / OfMIC]
  
# Remove unused columns to save space later
data.Ab <- data.Ab[, .(fct38, failure, totaldose_alloc, OfMIC, wt, sty1)]

# =====================
## Alternative model S1
# =====================
# outcome: FCT38
gam.fctS1.f <- bf(fct38 ~ t2(doseperMIC, DurRand, k=c(NULL, 4)) + s(wt) + s(sty1, bs="re"))
get_prior(
  gam.fctS1.f,
  data = data.Ab,
  family = "Gamma"
)
gam.Ab.fct.prior <- prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") +
  prior("gamma(2, 2)", class = "shape")

# get complete cases
data.Ab.fct <- data.Ab[!is.na(fct38)][!is.na(doseperMIC)][!is.na(DurRand)][!is.na(wt)][!is.na(sty1)]

# fit model
gam.Ab.fct <- brm(
  gam.fctS1.f,
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

## plot g-computation delineated by DurRand
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.fct, data.Ab.fct,
  xvar = "doseperMIC",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_fct_bydur_brms_S1.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$doseperMIC, na.rm = TRUE),
  xlab = "dose / MIC (mg/kg/day / mg/L)",
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## 2. outcome: failure
gam.failureS1.f <- bf(failure ~ t2(doseperMIC, DurRand, k=c(NULL, 4)) + s(wt) + s(sty1, bs="re"))
get_prior(
  gam.failureS1.f,
  data = data.Ab,
  family = "binomial"
)
gam.Ab.failure.prior <-prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") 
# get complete cases
data.Ab.failure <- data.Ab[!is.na(failure)][!is.na(doseperMIC)][!is.na(DurRand)][
  !is.na(wt)][!is.na(sty1)]
data.Ab.failure$failure <- as.numeric(data.Ab.failure$failure)

# fit model
gam.Ab.failure <- brm(
  gam.failureS1.f,
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
save(gam.Ab.fct, gam.Ab.failure, file = "./data/baseline_brms_models_S1.RData")
# load(file = "./data/baseline_brms_models_S1.RData")

## check posterior
pp_check(gam.Ab.failure)

## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.failure, data.Ab.failure,
  xvar = "doseperMIC",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_failure_bydur_brms_S1.pdf", gcomp_list,
  Edose = mean(data.Ab.failure$doseperMIC, na.rm = TRUE),
  ylab = "Expected change in\nfailure probability",
  xlab = "dose / MIC (mg/kg/day / mg/L)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

# =====================
## Alternative model S2
# =====================
# outcome: FCT38
gam.fctS2.f <- bf(fct38 ~ t2(dose, DurRand, OfMIC, k=c(NULL, 4, NULL)) + s(wt) + s(sty1, bs="re"))
get_prior(
  gam.fctS2.f,
  data = data.Ab,
  family = "Gamma"
)
gam.Ab.fct.prior <- prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") +
  prior("gamma(2, 2)", class = "shape")

# get complete cases
data.Ab.fct <- data.Ab[!is.na(fct38)][!is.na(dose)][!is.na(DurRand)][
  !is.na(OfMIC)][!is.na(wt)][!is.na(sty1)]

# fit model
gam.Ab.fct <- brm(
  gam.fctS2.f,
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

## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.fct, data.Ab.fct,
  xvar = "dose",
  bymic = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_fct_bymic_brms_S2.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  xlab = "dose (mg/kg/day)",
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## plot g-computation delineated by DurRand
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.fct, data.Ab.fct,
  xvar = "dose",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_fct_bydur_brms_S2b.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  xlab = "dose (mg/kg/day)",
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## 2. outcome: failure
gam.failureS2.f <- bf(
  failure ~ t2(dose, DurRand, OfMIC, k=c(NULL, 4, NULL)) + s(wt) + s(sty1, bs="re")
)
get_prior(
  gam.failureS2.f,
  data = data.Ab,
  family = "binomial"
)
gam.Ab.failure.prior <-prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") 
# get complete cases
data.Ab.failure <- data.Ab[!is.na(failure)][!is.na(dose)][!is.na(DurRand)][
  !is.na(OfMIC)][!is.na(wt)][!is.na(sty1)]
data.Ab.failure$failure <- as.numeric(data.Ab.failure$failure)

# fit model
gam.Ab.failure <- brm(
  gam.failureS2.f,
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
save(gam.Ab.fct, gam.Ab.failure, file = "./data/baseline_brms_models_S2.RData")

## check posterior
pp_check(gam.Ab.failure)

## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.failure, data.Ab.failure,
  xvar = "dose",
  bymic = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_failure_bymic_brms_S2.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  xlab = "dose (mg/kg/day)",
  ylab = "Expected change in\nfailure probability",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## plot g-computation delineated by Der
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.failure, data.Ab.failure,
  xvar = "dose",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_failure_bydur_brms_S2b.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  xlab = "dose (mg/kg/day)",
  ylab = "Expected change in\nfailure probability",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

# =====================
## Alternative model S3
# =====================
# outcome: FCT38
gam.fctS3.f <- bf(fct38 ~ t2(dose, DurRand, wt, k=c(NULL, 4, NULL)) + s(sty1, bs="re"))
get_prior(
  gam.fctS3.f,
  data = data.Ab,
  family = "Gamma"
)
gam.Ab.fct.prior <- prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") +
  prior("gamma(2, 2)", class = "shape")

# get complete cases
data.Ab.fct <- data.Ab[!is.na(fct38)][!is.na(dose)][!is.na(DurRand)][!is.na(wt)][!is.na(sty1)][!is.na(OfMIC)]

# fit model
gam.Ab.fct <- brm(
  gam.fctS3.f,
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

## plot g-computation delineated by DurRand
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.fct, data.Ab.fct,
  xvar = "dose",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_fct_bydur_brms_S3.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  xlab = "dose (mg/kg/day)",
  ylab = "Expected change in mean\nfever clearance time (h)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

## 2. outcome: failure
gam.failureS3.f <- bf(failure ~ t2(dose, DurRand, wt, k=c(NULL, 4, NULL)) + s(sty1, bs="re"))
get_prior(
  gam.failureS3.f,
  data = data.Ab,
  family = "binomial"
)
gam.Ab.failure.prior <-prior("normal(0, 0.5)", class = "b") + 
  prior("normal(0, 0.5)", class = "sds") +
  prior("normal(5, 1)", class = "Intercept") 
# get complete cases
data.Ab.failure <- data.Ab[!is.na(failure)][!is.na(dose)][!is.na(DurRand)][!is.na(wt)][!is.na(sty1)][!is.na(OfMIC)]
data.Ab.failure$failure <- as.numeric(data.Ab.failure$failure)

# fit model
gam.Ab.failure <- brm(
  gam.failureS3.f,
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
save(gam.Ab.fct, gam.Ab.failure, file = "./data/baseline_brms_models_S3.RData")
# load(file = "./data/baseline_brms_models.RData")

## check posterior
pp_check(gam.Ab.failure)

## plot g-computation delineated by MIC
gcomp_list <- calculate_gcomp_brms(
  gam.Ab.failure, data.Ab.failure,
  xvar = "dose",
  bydur = TRUE,
  refresh = 10
)
plot_gcomp_brms(
  "./figs/gformula_effect_failure_bydur_brms_S3.pdf", gcomp_list,
  Edose = mean(data.Ab.fct$dose, na.rm = TRUE),
  ylab = "Expected change in\nfailure probability",
  xlab = "dose (mg/kg/day)",
  width = 17 / 2.54, height = 17 / 3 / 2.54  
)

# =====================
## Compare everything
# =====================
load(file = "./data/baseline_brms_models_S3.RData")
looS3.failure <- loo(gam.Ab.failure)
looS3.fct <- loo(gam.Ab.fct)
load(file = "./data/baseline_brms_models_S2.RData")
looS2.failure <- loo(gam.Ab.failure)
looS2.fct <- loo(gam.Ab.fct)
load(file = "./data/baseline_brms_models_S1.RData")
looS1.failure <- loo(gam.Ab.failure)
looS1.fct <- loo(gam.Ab.fct)
load(file = "./data/baseline_brms_models.RData")
loo.failure <- loo(gam.Ab.failure)
loo.fct <- loo(gam.Ab.fct)

-2 * loo_compare(list(
  base = loo.fct,
  S1 = looS1.fct,
  S2 = looS2.fct,
  S3 = looS3.fct
))
-2 * loo_compare(list(
  base = loo.failure,
  S1 = looS1.failure,
  S2 = looS2.failure,
  S3 = looS3.failure
))

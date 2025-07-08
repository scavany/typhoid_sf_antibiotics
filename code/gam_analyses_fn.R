plot_counterfactuals <- function(file_name, data, model, outcome,
                                 exposure = "dose",
                                 cov1 = "DurRand",
                                 cov2 = "wt") {
  if (cov1 == "OfMIC") {
    DurRand.vec <- 2^c(-6, -3, 0)
  } else {
    DurRand.vec <- quantile(data.Ab[[cov1]], c(0.25, 0.5, 0.75), na.rm = TRUE)
  }
  if (cov2 == "OfMIC") {
    wt.vec <- 2^c(-6, -3, 0)
  } else {
    wt.vec <- quantile(data.Ab[[cov2]], c(0.25, 0.5, 0.75), na.rm = TRUE)
  }
  dose.vec <- seq(
    min(data.Ab[[exposure]], na.rm = TRUE),
    max(data.Ab[[exposure]], na.rm = TRUE), 
    0.1
  )
  pdf(
    file_name,
    width = 10,
    height = 5
  )
  par(mfrow=c(3,3),mar = c(4.1,4.1,2.1,2.1))
  if (outcome == "fct" || outcome == "fct38") {
    outcome.text <- "fever clearance time (h)"
  } else if (outcome == "failure") {
    outcome.text <- "failure probability"
  } else {
    outcome.text <- "OUTCOME"
  }
  if (exposure == "dose") {
    exposure.text <- "dose (mg/kg/d)"
  } else if (exposure == "doseperMIC") {
    exposure.text <- "dose / MIC (mg/kg/day / mg/L)"
  } else if (exposure == "totaldose_alloc") {
    exposure.text <- "total dose (mg/kg)"
  } else {
    exposure.text <- "EXPOSURE"
  }
  for (ii in seq_along(wt.vec)) {
    for (jj in seq_along(DurRand.vec)) {
      nd <- data.frame(
        dose.vec,
        rep(DurRand.vec[jj], length(dose.vec)),
        rep(wt.vec[ii], length(dose.vec)),
        rep(median(data.Ab$OfMIC, na.rm=TRUE), length(dose.vec)),
        rep(mean(data.Ab$Age.scaled, na.rm=TRUE), length(dose.vec)),
        rep(mean(data.Ab$Dur, na.rm=TRUE), length(dose.vec)),
        rep(mean(data.Ab$wt, na.rm=TRUE), length(dose.vec))
      )
      names(nd) <- c(exposure, cov1, cov2, "OfMIC", "Age.scaled", "Dur", "wt")
      prediction <- predict(
        model,
        newdata = nd,
        exclude = "s(sty1)",
        newdata.guaranteed = TRUE,
        type = "link",
        se.fit = TRUE
      )
      plot(
        dose.vec, inv_link(model)(prediction$fit),
        type = "l", lwd = 3,
        bty = "n",
        ylim = c(0, max(data.Ab[[outcome]], na.rm = TRUE)),
        xlab = ifelse(ii == 3, exposure.text, ""),
        ylab = ifelse(jj == 1, outcome.text, ""),
        main = paste0(
          cov1,
          " = ", DurRand.vec[jj],ifelse(cov1=="DurRand", "d", ""), ", ",
          cov2, " = ", wt.vec[ii], ifelse(cov2=="wt", "kg", "")
        )
      )
      abline(v = min(data.Ab[[exposure]], na.rm = TRUE), lty = 3)
      ## abline(h = fct.pt, lty = 3, col = "red")
      if (outcome == "fct38" || outcome == "fct") {
        ## if ((ii == 2) && (jj == 2)) {
        ##   text(
        ##     0.6 * 1000, 0.75* max(data.Ab[[outcome]], na.rm = TRUE),
        ##     paste0("Pre-treatment FCT: ",fct.pt,"h"),
        ##     col = "red", cex = 1.5
        ##   )
        ## }
        polygon(
          c(dose.vec,rev(dose.vec)),
          c(inv_link(model)(prediction$fit + 2 * prediction$se.fit),
            rev(inv_link(model)(pmax(1e-10,prediction$fit - 2 * prediction$se.fit)))),
          col=adjustcolor("black", 0.2), border=FALSE
        )
      } else {
        polygon(
          c(dose.vec,rev(dose.vec)),
          c(inv_link(model)(prediction$fit + 2 * prediction$se.fit),
            rev(inv_link(model)(prediction$fit - 2 * prediction$se.fit))),
          col=adjustcolor("black", 0.2), border=FALSE
        )
      }
    }
  }
  dev.off()
}

generate_matrices <- function(model.fct,
                              model.failure,
                              median.vals,
                              Ab.meanlog,
                              Ab.sdlog,
                              papi.range,
                              dose.range,
                              var2.range,
                              var3.range = NULL) {
  dose.name <- gsub(".range","",deparse(substitute(dose.range)))
  var2.name <- gsub(".range","",deparse(substitute(var2.range)))
  ## print(dose.name)
  ## print(var2.name)
  if (is.null(var3.range)) {
    fct.mat <- matrix(NA, nrow = length(dose.range), ncol = length(var2.range))
    fct.baseline.mat <- failure.baseline.mat <- failure.mat <- fct.mat
    for (ii in seq_along(dose.range)) {
      ## print(ii)
      for (jj in seq_along(var2.range)) {
        dose.baseline <- dose.range[ii] 
        dose.received <- dose.baseline * papi.range / 100
        nd <- data.frame(
          tempcol = rep(NA, length(dose.received))
        )
        for (val.name in names(median.vals)) {
          nd[[val.name]] <- median.vals[[val.name]]
        }
        nd[["tempcol"]] <- NULL
        nd[[dose.name]] <- dose.received
        nd[[var2.name]] <- var2.range[jj]
        prediction.fct <- inv_link(model.fct)(
          predict(
            model.fct,
            newdata = nd,
            type = "link",
            exclude = "s(sty1)",
            newdata.guaranteed = TRUE
          )
        )
        prediction.failure <- inv_link(model.failure)(
          predict(
            model.failure,
            newdata = nd,
            type = "link",
            exclude = "s(sty1)",
            newdata.guaranteed = TRUE
          )
        )
        weightedpred.fct <- prediction.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.failure <- prediction.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.fct.baseline <- prediction.fct * dnorm(papi.range, 100, 3)
        weightedpred.failure.baseline <- prediction.failure * dnorm(papi.range, 100, 3)
        fct.mat[ii, jj] <- sum(weightedpred.fct)
        failure.mat[ii, jj] <- sum(weightedpred.failure)
        fct.baseline.mat[ii, jj] <- sum(weightedpred.fct.baseline)
        failure.baseline.mat[ii, jj] <- sum(weightedpred.failure.baseline)
      }
    }
    out.list <- list(fct.mat, failure.mat, fct.baseline.mat, failure.baseline.mat)
  } else {
    var3.name <- gsub(".range","",deparse(substitute(var3.range)))
    fct.arr <- array(NA, dim = c(length(dose.range), length(var2.range), length(var3.range)))
    fct.baseline.arr <- failure.baseline.arr <- failure.arr <- fct.arr
    for (ii in seq_along(dose.range)) {
      for (jj in seq_along(var2.range)) {
        for (kk in seq_along(var3.range)) {
          dose.baseline <- dose.range[ii] 
          dose.received <- dose.baseline * papi.range / 100
          nd <- data.frame(
            tempcol = rep(NA, length(dose.received))
          )
          for (val.name in names(median.vals)) {
            nd[[val.name]] <- median.vals[[val.name]]
          }
          nd[["tempcol"]] <- NULL
          nd[[dose.name]] <- dose.received
          nd[[var2.name]] <- var2.range[jj]
          nd[[var3.name]] <- var3.range[kk]
          prediction.fct <- inv_link(model.fct)(
            predict(
              model.fct,
              newdata = nd,
              type = "link",
              exclude = "s(sty1)",
              newdata.guaranteed = TRUE
            )
          )
          prediction.failure <- inv_link(model.failure)(
            predict(
              model.failure,
              newdata = nd,
              type = "link",
              exclude = "s(sty1)",
              newdata.guaranteed = TRUE
            )
          )
          weightedpred.fct <- prediction.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
          weightedpred.failure <- prediction.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
          weightedpred.fct.baseline <- prediction.fct * dnorm(papi.range, 100, 3)
          weightedpred.failure.baseline <- prediction.failure * dnorm(papi.range, 100, 3)
          fct.arr[ii, jj, kk] <- sum(weightedpred.fct)
          failure.arr[ii, jj, kk] <- sum(weightedpred.failure)
          fct.baseline.arr[ii, jj, kk] <- sum(weightedpred.fct.baseline)
          failure.baseline.arr[ii, jj, kk] <- sum(weightedpred.failure.baseline)
        }
      }
    }
    out.list <- list(fct.arr, failure.arr, fct.baseline.arr, failure.baseline.arr)
  }
  names(out.list) <- c("fct","failure","fct.baseline","failure.baseline")
  return(out.list)
}

generate_matrices_se <- function(model.fct,
                                 model.failure,
                                 median.vals,
                                 Ab.meanlog,
                                 Ab.sdlog,
                                 papi.range,
                                 dose.range,
                                 var2.range,
                                 var3.range = NULL) {
  dose.name <- gsub(".range","",deparse(substitute(dose.range)))
  var2.name <- gsub(".range","",deparse(substitute(var2.range)))
  ## print(dose.name)
  ## print(var2.name)
  if (is.null(var3.range)) {
    fct.mat <- matrix(NA, nrow = length(dose.range), ncol = length(var2.range))
    fct.baseline.mat <- failure.baseline.mat <- failure.mat <- fct.mat
    fct.baseline.matlo <- failure.baseline.matlo <- failure.matlo <- fct.matlo <- fct.mat
    fct.baseline.mathi <- failure.baseline.mathi <- failure.mathi <- fct.mathi <- fct.mat
    for (ii in seq_along(dose.range)) {
      ## print(ii)
      for (jj in seq_along(var2.range)) {
        dose.baseline <- dose.range[ii] 
        dose.received <- dose.baseline * papi.range / 100
        nd <- data.frame(
          tempcol = rep(NA, length(dose.received))
        )
        for (val.name in names(median.vals)) {
          nd[[val.name]] <- median.vals[[val.name]]
        }
        nd[["tempcol"]] <- NULL
        nd[[dose.name]] <- dose.received
        nd[[var2.name]] <- var2.range[jj]
        link.fct <- predict(
          model.fct,
          newdata = nd,
          type = "link",
          exclude = "s(sty1)",
          newdata.guaranteed = TRUE,
          se.fit = TRUE
        )
        link.failure <- predict(
          model.failure,
          newdata = nd,
          type = "link",
          exclude = "s(sty1)",
          newdata.guaranteed = TRUE,
          se.fit = TRUE
        )
        prediction.fct <- inv_link(model.fct)(link.fct$fit)
        predlo.fct <- inv_link(model.fct)(pmax(1e-10, link.fct$fit - 1.96 * link.fct$se.fit))
        predhi.fct <- inv_link(model.fct)(link.fct$fit + 1.96 * link.fct$se.fit)
        prediction.failure <- inv_link(model.failure)(link.failure$fit)
        predlo.failure <- inv_link(model.failure)(link.failure$fit - 1.96 * link.failure$se.fit)
        predhi.failure <- inv_link(model.failure)(link.failure$fit + 1.96 * link.failure$se.fit)
        weightedpred.fct <- prediction.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpredlo.fct <- predlo.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpredhi.fct <- predhi.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.failure <- prediction.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpredlo.failure <- predlo.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpredhi.failure <- predhi.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.fct.baseline <- prediction.fct * dnorm(papi.range, 100, 3)
        weightedpredlo.fct.baseline <- predlo.fct * dnorm(papi.range, 100, 3)
        weightedpredhi.fct.baseline <- predhi.fct * dnorm(papi.range, 100, 3)
        weightedpred.failure.baseline <- prediction.failure * dnorm(papi.range, 100, 3)
        weightedpredlo.failure.baseline <- predlo.failure * dnorm(papi.range, 100, 3)
        weightedpredhi.failure.baseline <- predhi.failure * dnorm(papi.range, 100, 3)
        fct.mat[ii, jj] <- sum(weightedpred.fct)
        fct.matlo[ii, jj] <- sum(weightedpredlo.fct)
        fct.mathi[ii, jj] <- sum(weightedpredhi.fct)
        failure.mat[ii, jj] <- sum(weightedpred.failure)
        failure.matlo[ii, jj] <- sum(weightedpredlo.failure)
        failure.mathi[ii, jj] <- sum(weightedpredhi.failure)
        fct.baseline.mat[ii, jj] <- sum(weightedpred.fct.baseline)
        fct.baseline.matlo[ii, jj] <- sum(weightedpredlo.fct.baseline)
        fct.baseline.mathi[ii, jj] <- sum(weightedpredhi.fct.baseline)
        failure.baseline.mat[ii, jj] <- sum(weightedpred.failure.baseline)
        failure.baseline.matlo[ii, jj] <- sum(weightedpredlo.failure.baseline)
        failure.baseline.mathi[ii, jj] <- sum(weightedpredhi.failure.baseline)
      }
    }
    out.list <- list(fct.mat, failure.mat, fct.baseline.mat, failure.baseline.mat,
                     fct.matlo, failure.matlo, fct.baseline.matlo, failure.baseline.matlo,
                     fct.mathi, failure.mathi, fct.baseline.mathi, failure.baseline.mathi)
  } ## else {
  ##   var3.name <- gsub(".range","",deparse(substitute(var3.range)))
  ##   fct.arr <- array(NA, dim = c(length(dose.range), length(var2.range), length(var3.range)))
  ##   fct.baseline.arr <- failure.baseline.arr <- failure.arr <- fct.arr
  ##   for (ii in seq_along(dose.range)) {
  ##     for (jj in seq_along(var2.range)) {
  ##       for (kk in seq_along(var3.range)) {
  ##         dose.baseline <- dose.range[ii] 
  ##         dose.received <- dose.baseline * papi.range / 100
  ##         nd <- data.frame(
  ##           tempcol = rep(NA, length(dose.received))
  ##         )
  ##         for (val.name in names(median.vals)) {
  ##           nd[[val.name]] <- median.vals[[val.name]]
  ##         }
  ##         nd[["tempcol"]] <- NULL
  ##         nd[[dose.name]] <- dose.received
  ##         nd[[var2.name]] <- var2.range[jj]
  ##         nd[[var3.name]] <- var3.range[kk]
  ##         prediction.fct <- inv_link(model.fct)(
  ##           predict(
  ##             model.fct,
  ##             newdata = nd,
  ##             type = "link",
  ##             exclude = "s(sty1)",
  ##             newdata.guaranteed = TRUE
  ##           )
  ##         )
  ##         prediction.failure <- inv_link(model.failure)(
  ##           predict(
  ##             model.failure,
  ##             newdata = nd,
  ##             type = "link",
  ##             exclude = "s(sty1)",
  ##             newdata.guaranteed = TRUE
  ##           )
  ##         )
  ##         weightedpred.fct <- prediction.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
  ##         weightedpred.failure <- prediction.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
  ##         weightedpred.fct.baseline <- prediction.fct * dnorm(papi.range, 100, 3)
  ##         weightedpred.failure.baseline <- prediction.failure * dnorm(papi.range, 100, 3)
  ##         fct.arr[ii, jj, kk] <- sum(weightedpred.fct)
  ##         failure.arr[ii, jj, kk] <- sum(weightedpred.failure)
  ##         fct.baseline.arr[ii, jj, kk] <- sum(weightedpred.fct.baseline)
  ##         failure.baseline.arr[ii, jj, kk] <- sum(weightedpred.failure.baseline)
  ##       }
  ##     }
  ##   }
  ##   out.list <- list(fct.arr, failure.arr, fct.baseline.arr, failure.baseline.arr)
  ## }
  names(out.list) <- c("fct","failure","fct.baseline","failure.baseline",
                       "fct.low","failure.low","fct.low.baseline","failure.low.baseline",
                       "fct.high","failure.high","fct.high.baseline","failure.high.baseline")
  return(out.list)
}

generate_matrices_se_bootstrap <- function(data.Ab,
                                           model.fct,
                                           model.failure,
                                           median.vals,
                                           Ab.meanlog,
                                           Ab.sdlog,
                                           papi.range,
                                           dose.range,
                                           var2.range,
                                           R = 2^10,
                                           parallel = TRUE,
                                           type = "FORK",
                                           cores = 8,
                                           seed = 161) {
  if (parallel) {
    cl <- makeCluster(cores, type = type)
    registerDoParallel(cl)
    `%do.fn%` <- `%dopar%`
  } else {
    `%do.fn%` <- `%do%`
  }
  dose.name <- gsub(".range","",deparse(substitute(dose.range)))
  var2.name <- gsub(".range","",deparse(substitute(var2.range)))
  set.seed(seed)
  matlist <- foreach(icount(R), .packages = c("mgcv", "dplyr")) %do.fn% {
    fct.mat <- matrix(NA, nrow = length(dose.range), ncol = length(var2.range))
    fct.baseline.mat <- failure.baseline.mat <- failure.mat <- fct.mat
    indices <- sample(nrow(data.Ab), replace = TRUE)
    boot.sample <- data.Ab[indices, ]
    fit.boot.fct <- gam(
      formula(model.fct),
      data = boot.sample,
      family = family(model.fct),
      method = "REML"
    )
    fit.boot.failure <- gam(
      formula(model.failure),
      data = boot.sample,
      family = family(model.failure),
      method = "REML"
    )
    for (ii in seq_along(dose.range)) {
      for (jj in seq_along(var2.range)) {
        dose.baseline <- dose.range[ii]
        dose.received <- dose.baseline * papi.range / 100
        nd <- data.frame(
          tempcol = rep(NA, length(dose.received))
        )
        for (val.name in names(median.vals)) {
          nd[[val.name]] <- median.vals[[val.name]]
        }
        nd[["tempcol"]] <- NULL
        nd[[dose.name]] <- dose.received
        nd[[var2.name]] <- var2.range[jj]
        link.fct <- predict(
          fit.boot.fct,
          newdata = nd,
          type = "link",
          exclude = "s(sty1)",
          newdata.guaranteed = TRUE,
          se.fit = TRUE
        )
        link.failure <- predict(
          fit.boot.failure,
          newdata = nd,
          type = "link",
          exclude = "s(sty1)",
          newdata.guaranteed = TRUE,
          se.fit = TRUE
        )
        prediction.fct <- inv_link(fit.boot.fct)(link.fct$fit)
        prediction.failure <- inv_link(fit.boot.failure)(link.failure$fit)
        weightedpred.fct <- prediction.fct * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.failure <- prediction.failure * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
        weightedpred.fct.baseline <- prediction.fct * dnorm(papi.range, 100, 3)
        weightedpred.failure.baseline <- prediction.failure * dnorm(papi.range, 100, 3)
        fct.mat[ii, jj] <- sum(weightedpred.fct)
        failure.mat[ii, jj] <- sum(weightedpred.failure)
        fct.baseline.mat[ii, jj] <- sum(weightedpred.fct.baseline)
        failure.baseline.mat[ii, jj] <- sum(weightedpred.failure.baseline)
      }
    }
    out.list <- list(fct.mat, failure.mat, fct.baseline.mat, failure.baseline.mat)
    names(out.list) <- c("fct","failure","fct.baseline","failure.baseline")
    out.list
  }
  if (parallel) stopCluster(cl)
  names(matlist) <- 1:R
  return(combine_output_arrays(matlist))
}

combine_output_arrays <- function(z) {
  z <- lapply(z, `[`, names(z[[1]]))
  z <- apply(do.call(rbind, z), 2, as.list)
  lapply(z, function(x) abind(x, along = 3))
}

generate_matrices_faeces <- function(model,
                                     median.vals,
                                     Ab.meanlog,
                                     Ab.sdlog,
                                     papi.range,
                                     dose.range,
                                     var2.range) {
  dose.name <- gsub(".range","",deparse(substitute(dose.range)))
  var2.name <- gsub(".range","",deparse(substitute(var2.range)))
  faeces.mat <- matrix(NA, nrow = length(dose.range), ncol = length(var2.range))
  faeces.baseline.mat <- faeces.mat
  for (ii in seq_along(dose.range)) {
    for (jj in seq_along(var2.range)) {
      dose.baseline <- dose.range[ii] 
      dose.received <- dose.baseline * papi.range / 100
      nd <- data.frame(
        tempcol = rep(NA, length(dose.received))
      )
      for (val.name in names(median.vals)) {
        nd[[val.name]] <- median.vals[[val.name]]
      }
      nd[["tempcol"]] <- NULL
      nd[[dose.name]] <- dose.received
      nd[[var2.name]] <- var2.range[jj]
      prediction.faeces <- inv_link(model)(
        predict(
          model,
          newdata = nd,
          type = "link",
          exclude = "s(sty1)",
          newdata.guaranteed = TRUE
        )
      )
      weightedpred.faeces <- prediction.faeces * dlnorm(papi.range, Ab.meanlog, Ab.sdlog)
      weightedpred.faeces.baseline <- prediction.faeces * dnorm(papi.range, 100, 3)
      faeces.mat[ii, jj] <- sum(weightedpred.faeces)
      faeces.baseline.mat[ii, jj] <- sum(weightedpred.faeces.baseline)
    }
  }
  out.list <- list(faeces.mat, faeces.baseline.mat)
  names(out.list) <- c("faeces","faeces.baseline")
  return(out.list)
}

plot_sfeffect <- function(file_stem,
                          dose.range,
                          var2.range,
                          out.list,
                          xlab = "dose (mg/kg)",
                          ylab = "duration (days)") {
  png(paste0(file_stem, "_fct.png"))
  ## par(mar = c(4.1, 4.1, 1.1, 1.1))
  outcome <- out.list[["fct"]] - out.list[["fct.baseline"]]
  zlim.fct <- range(outcome, finite = TRUE)
  nlevels.fct <- 20
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  ## col.fct[1] <- "blue"
  filled.contour(
    dose.range,
    var2.range,
    outcome,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct
  )
  dev.off()
  png(paste0(file_stem, "_failure.png"))
  outcome <- out.list[["failure"]] - out.list[["failure.baseline"]]
  zlim.failure <- range(outcome, finite = TRUE)
  nlevels.failure <- 20
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  ## col.failure[1] <- "blue"
  filled.contour(
    dose.range,
    var2.range,
    outcome,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure
  )
  dev.off()
}

plot_sfeffect_se <- function(file_stem,
                             dose.range,
                             var2.range,
                             out.list,
                             Edose = 52,
                             xlab = "dose (mg/kg)",
                             ylab = "duration (days)") {
  ## par(mar = c(4.1, 4.1, 1.1, 1.1))
  outcome <- out.list[["fct"]] - out.list[["fct.baseline"]]
  outcome.lo <- out.list[["fct.low"]] - out.list[["fct.low.baseline"]]
  outcome.lo[which(outcome.lo > 10)] <- 10
  outcome.lo[which(outcome.lo < -10)] <- -10
  outcome.hi <- out.list[["fct.high"]] - out.list[["fct.high.baseline"]]
  zlim.fct <- range(outcome, finite = TRUE)
  zlim.fct.lo <- range(outcome.lo, finite = TRUE)
  zlim.fct.hi <- range(outcome.hi, finite = TRUE)
  nlevels.fct <- 20
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  levels.fct.lo <- pretty(zlim.fct.lo, nlevels.fct)
  col.fct.lo <- color.palette.fn(levels.fct.lo)
  levels.fct.hi <- pretty(zlim.fct.hi, nlevels.fct)
  col.fct.hi <- color.palette.fn(levels.fct.hi)
  ## col.fct[1] <- "blue"
  png(
    "temp_panel_1.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.lo,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct.lo,
    main = "at lower 95% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_2.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct,
    main = "at mean",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_3.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.hi,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct.hi,
    main = "at upper 95% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    paste0(file_stem, "_fct.png"),
    res = 600, width = 17 * 600 / 2.54,
    height = 8.5 * 600 / 2.54 * (5 / 6)
  )
  layout(
    matrix(
      c(
        rep(c(rep(0, 5), rep(1, 10), rep(0, 5)), 12),
        rep(2, 24 * 20),
        rep(c(rep(0, 5), rep(3, 10), rep(0, 5)), 12)
      ),
      nrow = 20
    )
  )
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  im1 <- load.image("temp_panel_1.png")
  im2 <- load.image("temp_panel_2.png")
  im3 <- load.image("temp_panel_3.png")
  plot(im1, axes = FALSE)
  plot(im2, axes = FALSE)
  plot(im3, axes = FALSE)
  dev.off()
  outcome <- out.list[["failure"]] - out.list[["failure.baseline"]]
  outcome.lo <- out.list[["failure.low"]] - out.list[["failure.low.baseline"]]
  outcome.hi <- out.list[["failure.high"]] - out.list[["failure.high.baseline"]] 
  zlim.failure <- range(outcome, finite = TRUE)
  zlim.failure.lo <- range(outcome.lo, finite = TRUE)
  zlim.failure.hi <- range(outcome.hi, finite = TRUE)
  nlevels.failure <- 20
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  levels.failure.lo <- pretty(zlim.failure.lo, nlevels.failure)
  levels.failure.hi <- pretty(zlim.failure.hi, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  col.failure.lo <- color.palette.fn(levels.failure.lo)
  col.failure.hi <- color.palette.fn(levels.failure.hi)
  ## col.failure[1] <- "blue"
  png(
    "temp_panel_1.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.lo,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure.lo,
    main = "at lower 95% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_2.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure,
    main = "baseline model",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_3.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.hi,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure.hi,
    main = "at upper 95% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    paste0(file_stem, "_failure.png"),
    res = 600, width = 17 * 600 / 2.54,
    height = 8.5 * 600 / 2.54 * (5 / 6)
  )
  layout(
    matrix(
      c(
        rep(c(rep(0, 5), rep(1, 10), rep(0, 5)), 12),
        rep(2, 24 * 20),
        rep(c(rep(0, 5), rep(3, 10), rep(0, 5)), 12)
      ),
      nrow = 20
    )
  )
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  im1 <- load.image("temp_panel_1.png")
  im2 <- load.image("temp_panel_2.png")
  im3 <- load.image("temp_panel_3.png")
  plot(im1, axes = FALSE)
  plot(im2, axes = FALSE)
  plot(im3, axes = FALSE)
  dev.off()
  system("rm temp_panel*png")
}

plot_sfeffect_se_bootstrap <- function(file_stem,
                                       dose.range,
                                       var2.range,
                                       out.list,
                                       matlist,
                                       Edose = 52,
                                       xlab = "dose (mg/kg)",
                                       ylab = "duration (days)") {
  ## par(mar = c(4.1, 4.1, 1.1, 1.1))
  outcome.baseline <-   pmax(out.list[["fct"]], 0.0) - pmax(out.list[["fct.baseline"]], 0.0)
  outcome <- pmax(matlist[["fct"]], 0.0) - pmax(matlist[["fct.baseline"]], 0.0)
  outcome.lo <- apply(
    outcome,
    c(1, 2),
    function(x) quantile(x, 0.05)
  )
  outcome.hi <- apply(
    outcome,
    c(1, 2),
    function(x) quantile(x, 0.95)
  )
  zlim.fct <- range(outcome.baseline, finite = TRUE)
  zlim.fct.lo <- range(outcome.lo, finite = TRUE)
  zlim.fct.hi <- range(outcome.hi, finite = TRUE)
  nlevels.fct <- 20
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  levels.fct.lo <- pretty(zlim.fct.lo, nlevels.fct)
  col.fct.lo <- color.palette.fn(levels.fct.lo)
  levels.fct.hi <- pretty(zlim.fct.hi, nlevels.fct)
  col.fct.hi <- color.palette.fn(levels.fct.hi)
  ## col.fct[1] <- "blue"
  png(
    "temp_panel_1.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.lo,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct.lo,
    main = "lower 90% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_2.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.baseline,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct,
    main = "baseline model",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_3.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.hi,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct.hi,
    main = "upper 90% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    paste0(file_stem, "_fct.png"),
    res = 600, width = 17 * 600 / 2.54,
    height = 8.5 * 600 / 2.54 * (5 / 6)
  )
  layout(
    matrix(
      c(
        rep(c(rep(0, 5), rep(1, 10), rep(0, 5)), 12),
        rep(2, 24 * 20),
        rep(c(rep(0, 5), rep(3, 10), rep(0, 5)), 12)
      ),
      nrow = 20
    )
  )
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  im1 <- load.image("temp_panel_1.png")
  im2 <- load.image("temp_panel_2.png")
  im3 <- load.image("temp_panel_3.png")
  plot(im1, axes = FALSE)
  plot(im2, axes = FALSE)
  plot(im3, axes = FALSE)
  dev.off()
  outcome.baseline <- out.list[["failure"]] - out.list[["failure.baseline"]]
  outcome <- matlist[["failure"]] - matlist[["failure.baseline"]]
  outcome.lo <- apply(
    outcome,
    c(1, 2),
    function(x) quantile(x, 0.05)
  )
  outcome.hi <- apply(
    outcome,
    c(1, 2),
    function(x) quantile(x, 0.95)
  )
  zlim.failure <- range(outcome.baseline, finite = TRUE)
  zlim.failure.lo <- range(outcome.lo, finite = TRUE)
  zlim.failure.hi <- range(outcome.hi, finite = TRUE)
  nlevels.failure <- 20
  ## ## OLD
  ## levels.failure <- pretty(zlim.failure, nlevels.failure)
  ## levels.failure.lo <- pretty(zlim.failure.lo, nlevels.failure)
  ## levels.failure.hi <- pretty(zlim.failure.hi, nlevels.failure)
  ## col.failure <- color.palette.fn(levels.failure)
  ## col.failure.lo <- color.palette.fn(levels.failure.lo)
  ## col.failure.hi <- color.palette.fn(levels.failure.hi)
  ## NEW
  zlim.failure <- range(
    outcome.baseline, -outcome.baseline,
    outcome.lo, -outcome.lo,
    outcome.hi, -outcome.hi
  )
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  png(
    "temp_panel_1.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.lo,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure, #.lo,
    zlim = zlim.failure,
    main = "lower 90% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_2.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.baseline,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure, 
    zlim = zlim.failure,
    main = "baseline model",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    "temp_panel_3.png",
    res = 600, width = 8.5 * 600 / 2.54, height = (5 / 6) * 8.5 * 600 / 2.54,
    pointsize = 8
  )
  par(mar = c(4.1, 4.1, 3.1, 0.1))
  filled.contour(
    dose.range,
    var2.range,
    outcome.hi,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure, #.hi,
    zlim = zlim.failure,
    main = "upper 90% CI",
    plot.axes = {
      axis(1,at = seq(0, 200, 50));
      axis(2,at = seq(0, 1, 0.2));
      abline(v = Edose, lty = 3, lwd = 2)
    }
  )
  dev.off()
  png(
    paste0(file_stem, "_failure.png"),
    res = 600, width = 17 * 600 / 2.54,
    height = 8.5 * 600 / 2.54 * (5 / 6)
  )
  layout(
    matrix(
      c(
        rep(c(rep(0, 5), rep(1, 10), rep(0, 5)), 12),
        rep(2, 24 * 20),
        rep(c(rep(0, 5), rep(3, 10), rep(0, 5)), 12)
      ),
      nrow = 20
    )
  )
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  im1 <- load.image("temp_panel_1.png")
  im2 <- load.image("temp_panel_2.png")
  im3 <- load.image("temp_panel_3.png")
  plot(im1, axes = FALSE)
  plot(im2, axes = FALSE)
  plot(im3, axes = FALSE)
  dev.off()
  system("rm temp_panel*png")
}

plot_sfeffect_faeces <- function(file_stem,
                                 dose.range,
                                 var2.range,
                                 out.list,
                                 xlab = "dose (mg/kg)",
                                 ylab = "duration (days)") {
  png(paste0(file_stem, "_faeces.png"))
  ## par(mar = c(4.1, 4.1, 1.1, 1.1))
  outcome <- out.list[["faeces"]] - out.list[["faeces.baseline"]]
  zlim.faeces <- range(outcome, finite = TRUE)
  nlevels.faeces <- 20
  levels.faeces <- pretty(zlim.faeces, nlevels.faeces)
  col.faeces <- color.palette.fn(levels.faeces)
  ## col.faeces[1] <- "blue"
  filled.contour(
    dose.range,
    var2.range,
    outcome,
    xlab = xlab,
    ylab = ylab,
    key.title = title(main = "faeces\npos diff"),
    nlevels = nlevels.faeces,
    col = col.faeces
  )
  dev.off()
}

plot_sfeffect_threeway <- function(file_stem,
                                   median.vals,
                                   dose.range,
                                   var2.range,
                                   var3.range,
                                   out.list,
                                   xlab = "dose (mg/kg)",
                                   v2lab = "duration (days)",
                                   v3lab = "MIC (mg/L)") {
  dose.name <- gsub(".range","",deparse(substitute(dose.range)))
  var2.name <- gsub(".range","",deparse(substitute(var2.range)))
  var3.name <- gsub(".range","",deparse(substitute(var3.range)))
  dose.medindex <- which.min(abs(dose.range - median.vals[[dose.name]]))
  var2.medindex <- which.min(abs(var2.range - median.vals[[var2.name]]))
  var3.medindex <- which.min(abs(var3.range - median.vals[[var3.name]]))
  ## First do FCT
  outcome <- out.list[["fct"]] - out.list[["fct.baseline"]]
  nlevels.fct <- 20
  png(paste0(file_stem, "_1_fct.png"))
  zlim.fct <- range(outcome[,,var3.medindex], finite = TRUE)
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  filled.contour(
    dose.range,
    var2.range,
    outcome[,,var3.medindex],
    xlab = xlab,
    ylab = v2lab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct
  )
  dev.off()
  png(paste0(file_stem, "_2_fct.png"))
  zlim.fct <- range(outcome[,var2.medindex,], finite = TRUE)
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  filled.contour(
    dose.range,
    var3.range,
    outcome[,var2.medindex,],
    xlab = xlab,
    ylab = v3lab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct
  )
  dev.off()
  png(paste0(file_stem, "_3_fct.png"))
  zlim.fct <- range(outcome[dose.medindex,,], finite = TRUE)
  levels.fct <- pretty(zlim.fct, nlevels.fct)
  col.fct <- color.palette.fn(levels.fct)
  filled.contour(
    var2.range,
    var3.range,
    outcome[dose.medindex,,],
    xlab = v2lab,
    ylab = v3lab,
    key.title = title(main = "fct\ndiff"),
    nlevels = nlevels.fct,
    col = col.fct
  )
  dev.off()
  ## Now failure
  outcome <- out.list[["failure"]] - out.list[["failure.baseline"]]
  nlevels.failure <- 20
  png(paste0(file_stem, "_1_failure.png"))
  zlim.failure <- range(outcome[,,var3.medindex], finite = TRUE)
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  filled.contour(
    dose.range,
    var2.range,
    outcome[,,var3.medindex],
    xlab = xlab,
    ylab = v2lab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure
  )
  dev.off()
  png(paste0(file_stem, "_2_failure.png"))
  zlim.failure <- range(outcome[,var2.medindex,], finite = TRUE)
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  filled.contour(
    dose.range,
    var3.range,
    outcome[,var2.medindex,],
    xlab = xlab,
    ylab = v3lab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure
  )
  dev.off()
  png(paste0(file_stem, "_3_failure.png"))
  zlim.failure <- range(outcome[dose.medindex,,], finite = TRUE)
  levels.failure <- pretty(zlim.failure, nlevels.failure)
  col.failure <- color.palette.fn(levels.failure)
  filled.contour(
    var2.range,
    var3.range,
    outcome[dose.medindex,,],
    xlab = v2lab,
    ylab = v3lab,
    key.title = title(main = "failure\nrisk diff"),
    nlevels = nlevels.failure,
    col = col.failure
  )
  dev.off()
}

color.palette.fn <- function(ls) {
  col.out <- rep(NA, length(ls) - 1)
  col.out[ls[-1] > 0] <- hcl.colors(sum(ls[-1] > 0), "Reds", rev = TRUE)
  col.out[!(ls[-1] > 0)] <- hcl.colors(sum(!(ls[-1] > 0)), "Blues")
  return(col.out)
}

gcomp.bootstrap <- function(model, data, R, dose.vec, Edose.index = 1, dose.name = "dose", MIC = NA,
                            parallel = TRUE, type = "FORK", cores = 8, seed = 161) {
  if (parallel) {
    cl <- makeCluster(cores, type = type)
    registerDoParallel(cl)
    `%do.fn%` <- `%dopar%`
  } else {
    `%do.fn%` <- `%do%`
  }
  set.seed(seed)
  pred.samples <- foreach(
    icount(R),
    .combine = rbind,
    .packages = c("mgcv", "dplyr")
  ) %do.fn% {
    pred.samples.temp <- rep(0, length(dose.vec))
    indices <- sample(nrow(data), replace = TRUE)
    boot.sample <- data[indices, ]
    fit.boot <- gam(
      formula(model),
      data = boot.sample,
      family = family(model),
      method = "REML"
    )
    boot.sample.new <- select(boot.sample, !all_of(dose.name))
    if (!is.na(MIC)) boot.sample.new$OfMIC <- MIC
    boot.sample.new[[dose.name]] <- dose.vec[Edose.index]
    pred.base <- predict(
      fit.boot, 
      newdata = boot.sample.new, 
      type = "response"
    )
    for (jj in seq_len(length(dose.vec))[-Edose.index]) {
      boot.sample.new[[dose.name]] <- dose.vec[jj]
      pred.new <- predict(
        fit.boot, 
        newdata = boot.sample.new, 
        type = "response"
      )
      pred.samples.temp[jj] <- mean(pred.new - pred.base, na.rm = TRUE)
    }
    pred.samples.temp
  }
  if (parallel) stopCluster(cl)
  rownames(pred.samples) <- NULL
  return(pred.samples)
}

gcomp.bootstrap.sequential <- function(model, data, R, dose.vec, dose.name = "dose",
                                       Edose.index = 1, seed = 161) {
  set.seed(seed)
  pred.samples <- matrix(0, nrow = R, ncol = length(dose.vec))
  for(ii in 1:R) {
    indices <- sample(nrow(data), replace = TRUE)
    boot.sample <- data[indices, ]
    boot.sample.new <- select(boot.sample, !all_of(dose.name))
    fit.boot <- gam(
      formula(model),
      data = boot.sample,
      family = family(model),
      method = "REML"
    )
    boot.sample.new[[dose.name]] <- dose.vec[Edose.index]
    pred.base <- predict(
      fit.boot, 
      newdata = boot.sample.new, 
      type = "response"
    )
    for (jj in seq_len(length(dose.vec))[-Edose.index]) {
      boot.sample.new[[dose.name]] <- dose.vec[jj]
      pred.new <- predict(
        fit.boot, 
        newdata = boot.sample.new, 
        type = "response"
      )
      pred.samples[ii, jj] <- mean(pred.new - pred.base, na.rm = TRUE)
    }
  }
  return(pred.samples)
}

create_posterior_pred_faeces <- function(data, draws) {
  Xs <- as.matrix(data$Xs) # ndata x neffects
  bscols <- which(grepl("bs", names(draws))) 
  bs <- as.matrix(draws[, bscols]) # nsamples x neffects
  # Zs_1_1
  Zs_1_1 <- as.matrix(data$Zs_1_1) # ndata x nknots
  sds_1cols <- which(grepl("sds_1", names(draws)) & grepl("[1]", names(draws))) 
  sds_1 <- unlist(draws[, sds_1cols[1]]) # nsamples
  zs_1_1cols <- which(grepl("zs_1_1", names(draws)))
  zs_1_1 <- as.matrix(draws[, zs_1_1cols]) # nsamples x nknots
  s_1_1 <- sds_1 * zs_1_1 # nsamples x nknots
  # Zs_1_2
  Zs_1_2 <- as.matrix(data$Zs_1_2) # ndata x nknots
  sds_1cols <- which(grepl("sds_1", names(draws)) & grepl("[1]", names(draws))) 
  sds_1 <- unlist(draws[, sds_1cols[2]]) # nsamples
  zs_1_2cols <- which(grepl("zs_1_2", names(draws)))
  zs_1_2 <- as.matrix(draws[, zs_1_2cols]) # nsamples x nknots
  s_1_2 <- sds_1 * zs_1_2 # nsamples x nknots
  # Zs_1_3
  Zs_1_3 <- as.matrix(data$Zs_1_3) # ndata x nknots
  sds_1cols <- which(grepl("sds_1", names(draws)) & grepl("[1]", names(draws))) 
  sds_1 <- unlist(draws[, sds_1cols[3]]) # nsamples
  zs_1_3cols <- which(grepl("zs_1_3", names(draws)))
  zs_1_3 <- as.matrix(draws[, zs_1_3cols]) # nsamples x nknots
  s_1_3 <- sds_1 * zs_1_3 # nsamples x nknots
  # Zs_2_1
  Zs_2_1 <- as.matrix(data$Zs_2_1) # ndata x nknots
  sds_2cols <- which(grepl("sds_2", names(draws)) & grepl("[1]", names(draws))) 
  sds_2 <- unlist(draws[, sds_2cols[1]]) # nsamples
  zs_2_1cols <- which(grepl("zs_2_1", names(draws)))
  zs_2_1 <- as.matrix(draws[, zs_2_1cols]) # nsamples x nknots
  s_2_1 <- sds_2 * zs_2_1 # nsamples x nknots
  # Zs_3_1
  Zs_3_1 <- as.matrix(data$Zs_3_1) # ndata x nknots
  sds_3cols <- which(grepl("sds_3", names(draws)) & grepl("[1]", names(draws))) 
  sds_3 <- unlist(draws[, sds_3cols[1]]) # nsamples
  zs_3_1cols <- which(grepl("zs_3_1", names(draws)))
  zs_3_1 <- as.matrix(draws[, zs_3_1cols]) # nsamples x nknots
  s_3_1 <- sds_3 * zs_3_1 # nsamples x nknots
  # all three outputs below have dim nsamples x ndata
  mu <- draws$Intercept + t(Xs %*% t(bs)) + t(Zs_1_1 %*% t(s_1_1)) + t(Zs_1_2 %*% t(s_1_2)) +
    t(Zs_1_3 %*% t(s_1_3)) + t(Zs_2_1 %*% t(s_2_1)) + t(Zs_3_1 %*% t(s_3_1)) 
  prob <- plogis(mu)
  pos.pred <- t(matrix(rbinom(length(prob), data$trials, t(prob)), ncol = nrow(draws)))
  return(list(mu = mu, p = prob, post = pos.pred))
}

calculate_gcomp_brms <- function(model, data,
                                 quantiles_desired = c(0.05, 0.95),
                                 xvar = "totaldose_alloc",
                                 bymic = FALSE,
                                 bydur = FALSE,
                                 byperc = FALSE,
                                 refresh = 0) {
  if (bymic & bydur) {
    stop("choose either bymic or bydur or neither; not both")
  }
  if (byperc) {
    dose.vec <- seq(0, 200, 5)
  } else {
    dose.vec <- seq(
      floor(min(data[[xvar]])),
      ceiling(max(data[[xvar]])),
      length.out = 150
    )
  }
  if (bymic) {
    MIC.vec <- 2^c(-3, -2, 0)# 2^c(-6, -3, 0)
  } else if (bydur) {
    MIC.vec <- c(2, 3, 5)
  } else {
    MIC.vec <- NULL
  }
  # set up data structures
  if (!bymic & !bydur) {
    mu_ate <- data.frame(matrix(
      NA, nrow = length(dose.vec), ncol = 1 + length(quantiles_desired)
    ))
    colnames(mu_ate) <- c("mean", quantiles_desired)
  } else {
    mu_ate <- array(
      NA,
      dim = c(length(MIC.vec), length(dose.vec), 1 + length(quantiles_desired)),
      dimnames = list(MIC.vec, dose.vec, c("mean", quantiles_desired))
    )
  }
  post_ate <- mu_ate
  # create baseline pred at mean xvar
  Edose <- mean(data[[xvar]], na.rm=TRUE)
  data.temp <- copy(data)
  if (!bymic & !bydur) {
    if (!byperc) data.temp[[xvar]] <- Edose
    mu_base <- posterior_epred(model, newdata = data.temp)
    post_base <- posterior_predict(model, newdata = data.temp)
    # now try for the range of counterfactual doses
    for (ii in seq_along(dose.vec)) {
      if (refresh > 0 & ii %% refresh == 0) print(ii)
      if (byperc) data.temp[[xvar]] <- data[[xvar]] * dose.vec[ii] / 100
      else data.temp[[xvar]] <- dose.vec[ii] 
      mu_temp <- posterior_epred(model, newdata = data.temp)
      post_temp <- posterior_predict(model, newdata = data.temp)
      mu_ate.temp <- rowMeans(mu_temp - mu_base)
      mu_ate[ii, 1] <- mean(mu_ate.temp)
      mu_ate[ii, -1] <- quantile(mu_ate.temp, probs = quantiles_desired)
      post_ate.temp <- rowMeans(post_temp - post_base)
      post_ate[ii, 1] <- mean(post_ate.temp)
      post_ate[ii, -1] <- quantile(post_ate.temp, probs = quantiles_desired)
    }
  } else if (bymic) {
    for (jj in seq_along(MIC.vec)) {
      data.temp$OfMIC <- MIC.vec[jj]
      # create baseline pred at mean dose (53.7) and MIC of interest
      if (!byperc) data.temp[[xvar]] <- Edose
      mu_base <- posterior_epred(model, newdata = data.temp)
      post_base <- posterior_predict(model, newdata = data.temp)
      for (ii in seq_along(dose.vec)) {
        if (refresh > 0 && (ii %% refresh == 0)) print(ii)
        if (byperc) data.temp[[xvar]] <- data[[xvar]] * dose.vec[ii] / 100
        else data.temp[[xvar]] <- dose.vec[ii] 
        mu_temp <- posterior_epred(model, newdata = data.temp)
        post_temp <- posterior_predict(model, newdata = data.temp)
        mu_ate.temp <- rowMeans(mu_temp - mu_base)
        mu_ate[jj, ii, 1] <- mean(mu_ate.temp)
        mu_ate[jj, ii, -1] <- quantile(mu_ate.temp, probs = quantiles_desired)
        post_ate.temp <- rowMeans(post_temp - post_base)
        post_ate[jj, ii, 1] <- mean(post_ate.temp)
        post_ate[jj, ii, -1] <- quantile(post_ate.temp, probs = quantiles_desired)
      }
    }
  } else {
    for (jj in seq_along(MIC.vec)) {
      data.temp$DurRand <- MIC.vec[jj]
      # create baseline pred at mean dose (53.7) and MIC of interest
      if (!byperc) data.temp[[xvar]] <- Edose
      mu_base <- posterior_epred(model, newdata = data.temp)
      post_base <- posterior_predict(model, newdata = data.temp)
      for (ii in seq_along(dose.vec)) {
        if (refresh > 0 && (ii %% refresh == 0)) print(ii)
        if (byperc) data.temp[[xvar]] <- data[[xvar]] * dose.vec[ii] / 100
        else data.temp[[xvar]] <- dose.vec[ii] 
        mu_temp <- posterior_epred(model, newdata = data.temp)
        post_temp <- posterior_predict(model, newdata = data.temp)
        mu_ate.temp <- rowMeans(mu_temp - mu_base)
        mu_ate[jj, ii, 1] <- mean(mu_ate.temp)
        mu_ate[jj, ii, -1] <- quantile(mu_ate.temp, probs = quantiles_desired)
        post_ate.temp <- rowMeans(post_temp - post_base)
        post_ate[jj, ii, 1] <- mean(post_ate.temp)
        post_ate[jj, ii, -1] <- quantile(post_ate.temp, probs = quantiles_desired)
      }
    }
  }
  return(
    list(
      dose.vec = dose.vec, MIC.vec = MIC.vec,
      mu_ate = mu_ate, post_ate = post_ate
    )
  )
}

plot_gcomp_brms <- function(filename, gcomp_list,
                            Edose = 52,
                            xlab = "total dose (mg/kg)",
                            ylab = "Expected change in mean fever clearance time (h)",
                            width = 12 / 2.54, height = 12 / 2.54, pointsize = 13) {
  pdf(filename, width = width, height = height, pointsize = pointsize)
  if (is.null(gcomp_list$MIC.vec)) {  
    plot(
      gcomp_list$dose.vec,
      gcomp_list$mu_ate[, "mean"],
      ylim = range(gcomp_list$mu_ate[, c("0.05", "0.95")]),
      lty = 1,
      lwd = 2.5,
      type = "l",
      bty = "n",
      xaxs = "i",
      yaxs = "i",
      xlab = xlab,
      ylab = ylab
    )
    lines(gcomp_list$dose.vec, gcomp_list$mu_ate[, "0.05"], lty = 3)
    lines(gcomp_list$dose.vec, gcomp_list$mu_ate[, "0.95"], lty = 3)
    polygon(
      c(gcomp_list$dose.vec, rev(gcomp_list$dose.vec)),
      c(gcomp_list$mu_ate[, "0.05"], rev(gcomp_list$mu_ate[, "0.95"])),
      col = adjustcolor("grey",0.3), border=FALSE
    )
    abline(h = 0, lty = 2)
    abline(v = Edose, lty = 2)
    dev.off() 
  } else {
    par(mfrow = c(1, 3), mar = c(4.1, 5.1, 2.1, 1.1))
    for (ii in seq_along(gcomp_list$MIC.vec)) {
      if (min(gcomp_list$MIC.vec) < 1) {
        main.text <- ifelse(
          gcomp_list$MIC.vec[ii] == 1,
          "MIC: 1 mg/L",
          paste0("MIC: ", round(gcomp_list$MIC.vec[ii], 4), " mg/L")
        )
      } else {
        main.text <- paste0("Treatment duration: ", as.integer(gcomp_list$MIC.vec[ii]), " h")
      }
      plot(
        gcomp_list$dose.vec,
        gcomp_list$mu_ate[ii, , "mean"],
        ylim = range(gcomp_list$mu_ate[, , c("0.05", "0.95")]),
        lty = 1,
        lwd = 2.5,
        type = "l",
        bty = "n",
        xaxs = "i",
        yaxs = "i",
        main = main.text,
        xlab = xlab,
        ylab = ifelse(ii == 1, ylab, "")
      )
      lines(gcomp_list$dose.vec, gcomp_list$mu_ate[ii, , "0.05"], lty = 3)
      lines(gcomp_list$dose.vec, gcomp_list$mu_ate[ii, , "0.95"], lty = 3)
      polygon(
        c(gcomp_list$dose.vec, rev(gcomp_list$dose.vec)),
        c(gcomp_list$mu_ate[ii, , "0.05"], rev(gcomp_list$mu_ate[ii, , "0.95"])),
        col = adjustcolor("grey",0.3), border=FALSE
      )
      abline(h = 0, lty = 2)
      abline(v = Edose, lty = 2)
    }
    dev.off()
  }
}

plot_gcomp_brms_escmid <- function(filename, gcomp_list,
                                   Edose = 52,
                                   xlab = "total dose (mg/kg)",
                                   ylab = "Expected change in mean fever clearance time (h)",
                                   width = 12 / 2.54, height = 12 / 2.54, pointsize = 13) {
  pdf(filename, width = width, height = height, pointsize = pointsize)
  par(mfrow = c(1, 2), mar = 0.1 + c(4, 5, 2, 1))
  for (ii in 2:3) {
    if (min(gcomp_list$MIC.vec) < 1) {
      main.text <- "" # ifelse(
        # gcomp_list$MIC.vec[ii] == 1,
        # "Resistant\nMIC: 1 mg/L",
        # paste0("Sensitive\nMIC: 1/", as.integer(1/gcomp_list$MIC.vec[ii]), " mg/L")
      # )
    } else {
      main.text <- "" # paste0("Treatment duration: ", as.integer(gcomp_list$MIC.vec[ii]), " h")
    }
    plot(
      gcomp_list$dose.vec,
      gcomp_list$mu_ate[ii, , "mean"],
      ylim = range(gcomp_list$mu_ate[2:3, , c("0.05", "0.95")]),
      lty = 1,
      lwd = 5,
      type = "l",
      bty = "n",
      xaxs = "i",
      yaxs = "i",
      main = main.text,
      xlab = xlab,
      ylab = ifelse(ii == 2, ylab, "")
    )
    # mtext(ifelse(ii == 2, ylab, ""), 2, 3)
    lines(gcomp_list$dose.vec, gcomp_list$mu_ate[ii, , "0.05"], lty = 3, lwd = 2.5)
    lines(gcomp_list$dose.vec, gcomp_list$mu_ate[ii, , "0.95"], lty = 3, lwd = 2.5)
    polygon(
      c(gcomp_list$dose.vec, rev(gcomp_list$dose.vec)),
      c(gcomp_list$mu_ate[ii, , "0.05"], rev(gcomp_list$mu_ate[ii, , "0.95"])),
      col = adjustcolor("grey",0.3), border=FALSE
    )
    abline(h = 0, lty = 2)
    abline(v = Edose, lty = 2)
  }
  dev.off()
}

### Functions for the manuscript:
# Optimising time samples for determining AUC of pharmacokinetic data using 
# non-compartmental analysis
# Hughes JH, Upton RN, Reuter SE, Phelps MA, Foster DJR
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pred.sumexp <- function(par, x, d = 0) {
# Provides predictions according to model parameters
# par = sum of exponential parameters
# x = independent variable (time)
# d = order of derivative (uses dth derivative of model)
# Define objects
  l <- length(par)  # number of parameters
  a <- l %% 2 == 1  # absorption status (odd parameter length == absorption)
  n <- ceiling(l/2)  # number of exponentials
  m <- -abs(par[1:n])  # slope parameters (prevents exponential growth)
  b <- par[(n+1):l]  # intercept parameters
# Order parameters (allows for flip-flop)
  m.ord <- order(m, decreasing = T)  # slope order (terminal first, absorption last)
  b.ord <- m.ord  # intercept order (match slopes)
  if (a) b.ord <- order(m[-m.ord[n]], decreasing = T)  # if absorption curve remove extra term
  p <- c(m[m.ord], b[b.ord])  # ordered parameters
# Sum of exponentials
  for (i in 1:n) {  # for each exponential
    if (i == 1) yhat <- p[i]^d*exp(p[i]*x + p[n+i])  # first exponential defines yhat
    else if (i != n | !a) yhat <- yhat + p[i]^d*exp(p[i]*x + p[n+i])  # following exponentials add to yhat
    else if (a) yhat <- yhat - p[i]^d*exp(p[i]*x)*sum(exp(p[(n+1):(2*n-1)]))  # for absorption curve apply final term
  }
  return(yhat)  # predicted dependent variable (drug concentration)
}

order.sumexp <- function(par, n, a) {
# Sorts sum of exponential parameters to enable pred.sumexp
# par = sum of exponential parameters
# n = number of exponentials
# a = absorption status
  m <- -abs(head(par, n+a))  # slope parameters (prevents exponential growth)
  b <- tail(par, -(n+a))  # intercept parameters
  m.ord <- order(m, decreasing = T)  # slope order (terminal first, absorption last)
  b.ord <- m.ord  # intercept order (match slopes)
  if (a) b.ord <- order(m[-m.ord[n+a]], decreasing = T)  # if absorption curve remove extra term
  unname(c(m[m.ord], b[b.ord]))  # ordered parameters
}

mle.sumexp <- function(par, x, y, sigma, ga = F) {
# Determine log likelihood of given model parameters
# This function is used within optim.sumexp and is not recommended to be used 
# outside of that context
# par = sum of exponential parameters
# x = independent variable (time)
# y = observed dependent variable (drug concentration)
# sigma = proportional error
# ga = genetic algorithm status
  z <- ifelse(ga, 2, -2)  # adjust ofv for minimisation or maximisation
  yhat <- pred.sumexp(par, x)  # sum of exponential model prediction
  loglik <- dnorm(y, yhat, abs(yhat)*sigma, log = T)  # log likelihood
  return(z*sum(loglik))  # objective function value
}

pred.lambdaz <- function(dv, t) {
# Determine the slope of the terminal phase
# dv = observed dependent variable (drug concentration)
# t = independent variable (time)
# set up objects used iterative loop
# i = number of concentrations to be used to determine terminal slope
# begins with 2 as a safety, only uses slope from 2 if there are no more terminal concs
  i <- 2
  j <- 0
  bestR2 <- -1
  bestk <- rep(0, 3)
# Define the terminal concentrations as those after the maximum concentration
  terminal <- which(dv == max(dv))[1]:length(dv)
  rem <- matrix(length(dv)+1)
  if (length(terminal) >= i) {  # if more than 1 terminal concentration
    repeat {  # begin iterations
      for (l in 1:ncol(rem)) {  # for all combinations of i concentrations
      # find terminal slope for i concentrations and record linear model parameters
        mod <- suppressWarnings(lm(log(tail(dv[-rem[,l]], i)) ~ tail(unique(t)[-rem[,l]], i)))
        k <- unname(mod$coefficients)
        R2 <- suppressWarnings(as.numeric(summary(mod)["adj.r.squared"]))
        if (!is.finite(k[2])) {  # if slope is not finite break loop
          R2 <- -1
          break
        }
      # if adjusted r2 is NaN, use r2 instead
        if (is.nan(R2)) R2 <- suppressWarnings(as.numeric(summary(mod)["r.squared"]))
        if (k[2] < 0) {  # if terminal slope is negative
          if (R2 > bestR2) {  # if r2 is better than the best r2
            if (i > 2) bestR2 <- R2  # if number of concentrations in slope > 2
            bestk <- c(k, sd(residuals(mod)))
          } else {
            break  # break out of for loop
          }
        }
      }
      if (R2 < bestR2) break  # take model that is better than latest model
      if (i == length(terminal)) {  # once all points have been included test for model validity
        if (!is.finite(bestk[2])) { browser() }
        else if (bestk[2] == 0) {  # if all current models are invalid, trial removing random points
          if (length(terminal) != j + 3) {  # as long as there are still points to remove
            j <- j + 1
            i <- j + 2
            rem <- length(dv) - combn(i + 1, j) + 1
          } else {  # if not make a last ditch effort (intended for simulation study only)
            bestk <- c(max(dv), -log(2)/56, 0.01)
            break  # break out of repeat loop
          }
        } else { break }  # if there was a valid model then finish
      }
      i <- i + 1
    }
  }
  bestk
}

optim.sumexp <- function(data, oral = F, nexp = 3) {
# Determines best sum of exponential for given data
# data = mean pooled data;
#        first column independent variable;
#        second column dependent variable
# oral = whether the drug displays an absorption curve
# nexp = maximum fitted number of exponentials
# Set up objects
  res <- list(par = NULL, sumexp = NULL, value = NULL,
  error = NULL, hessian = NULL, message = NULL)
# Prepare data (remove t == 0, remove negative concentrations)
  x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
  y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
# Estimate candidate model parameters
  init <- pred.lambdaz(y, x)  # determine terminal slope for initial estimates
  for (i in 1:nexp) {  # for each of the number of exponentials
  # use genetic algorithms to search for good starting values for optimisation
    gares <- try(ga("real-valued",
      mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
      min = c(rep(init[2]*50, i + oral), rep(init[1]-2, i)),
      max = c(rep(init[2]/50, i + oral), rep(init[1]+2, i)),
      selection = gareal_lrSelection,
      crossover = gareal_spCrossover,
      mutation = gareal_raMutation,
      maxiter = 50,
      popSize = 250,
      monitor = F
    ))
    if (class(gares) == "try-error") {
      stop("Unsuccessful genetic algorithm optimisation: may be a bad set of random numbers please run again.")
    }
  # take initial parameters from genetic algorithm and use them for quasi-newton optimisation
    ga.par <- gares@solution[1, ]
    optres <- try(optim(
      ga.par,
      mle.sumexp,
      method = "BFGS", hessian = T,
      x = x, y = y, sigma = 0.01
    ))
    if (class(optres) == "try-error") {  # if quasi-Newton fails, use genetic algorithm result
      optres <- list(
        par = gares@solution[1,],
        value = -gares@fitnessValue,
        counts = NA,
        hessian = matrix(NA, ncol = length(ga.par), nrow = length(ga.par)),
        convergence = NA,
        message = NA
      )
    }
  # Create output
    par.ord <- order.sumexp(optres$par, i, oral)
    res$par[[i]] <- optres$par
    res$sumexp[[i]] <- par.ord
    res$value[[i]] <- c(ofv = optres$value, optres$counts)
    res$error[[i]] <- c(0.01, "fixed")
    res$hessian[[i]] <- optres$hessian
    res$message[[i]] <- c(convergence = optres$convergence,
      message = ifelse(is.null(optres$message), "NULL", optres$message)
    )
  }
  res
}

best.sumexp.aic <- function(opt) {
# Determines best sum of exponential parameters from optim.sumexp output
# Uses Aikaike's Information Criterion to determine best set of parameters
# opt = optim.sumexp output
  values <- unlist(opt$value)
  ofv <- values[which(names(values) == "ofv")]
  k <- unlist(lapply(opt$par, length))
  aic <- ofv + 2*k
  return(sapply(opt, function(x) x[which(aic == min(aic))]))
}

tmax.sumexp <- function(par, tlast = 24, res = 0.1) {
# Determines tmax for a set of sum of exponential parameters
# par = sum of exponential parameters
# tlast = final time checked for tmax (if tmax = tlast, tlast should be increased)
# res = resolution of search for tmax (default is 0.1 time units)
  times <- seq(0, tlast, by = res)
  yhat <- pred.sumexp(par, times)
  return(times[which(yhat == max(yhat))])
}

err.interv <- function(par, exp.par, tfirst, tlast, a = F, tmax = NULL) {
# Determines the trapezoidal error associated with a set of times for a
# specific sum of exponential equation
# This function is used within optim.sumexp and is not recommended to be used 
# outside of that context
# par = sample times for which trapezoidal error is calculated
# exp.par = sum of exponential parameters
# tfirst = first time for time samples (usually 0) 
# tlast = final time for time samples
# a = logical statement of whether exp.par is for an absorption curve
# tmax = tmax to be fixed if using that alternate method
# First put together times from first and last times, tmax and time intervals
  times <- c(tfirst, cumsum(par), tmax, tlast)
  times <- sort(times)
# Calculated the change in time between times
  deltat <- diff(times)
# Determine highest second derivative value for each interval
  if (a) {
  # if an absorption curve then must determine whether the first or second value 
  # is associated with the most error
    all.secd <- abs(pred.sumexp(exp.par, times, 2))
    secd <- c(NULL)
    for (i in 1:(length(times)-1)) {
      secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
    }
  } else {  
  # otherwise highest second derivative is always for the first time point
    secd <- pred.sumexp(exp.par, times[-length(times)], 2)
  }
# Determine trapezoidal error
  err <- deltat^3*secd/12
  sumerr <- sqrt(mean(err^2))
  return(sumerr)
}

optim.interv <- function(par, times, tmax = FALSE) {
# Determines best set of intervals to an sum of exponential equation
# par = sum of exponential equation parameters (typically output of best.sumexp)
# times = times to be optimised, should be the same length as the output you want
# should also define the first and last time to be measured
# tmax = logical statement of whether the fixed tmax alternate method should be used
# Set up environment by defining often used objects
  tfirst <- min(times)
  tlast <- max(times)
  tbound <- tlast - tfirst
  absorp <- ifelse((length(par) %% 2) != 0, T, F)
  npar <- length(times) - (2 + tmax*absorp)
# Calculated tmax if method is fixed tmax method is used and sum of exponentials 
# describes an absorption curve
  if (tmax & absorp) tmax.val <- tmax.sumexp(par, tlast, tlast/720)  # maximum length of 721
  else tmax.val <- NULL
# Set up nrep object to count number of iterations
  nrep <- 0
  repeat {
  # create randomised initial parameters for optimisation
    init.par <- cumsum(rep(tbound/2^runif(1, 4, 6), npar))
    res <- try(optim(
      init.par,
      err.interv,
      method = "L-BFGS-B", hessian = T,
      lower = tbound/(0.96*npar^2), upper = tbound/(npar/1.5), exp.par = par,
      tfirst = tfirst, tlast = tlast, a = absorp, tmax = tmax.val
    ))
    if (class(res) == "try-error") {stop("Critical error: bad seed or bad input")}
  # if optimisation converges or has run 10 times break
    if (res$convergence == 0 | nrep == 10) break
    nrep <- nrep + 1
  }
  res$times <- sort(c(cumsum(c(tfirst, res$par)), tmax.val, tlast))
  return(res)
}

obs.tlast <- function(obs) {
# Determines the best tlast for a set of observed data by calculating the slope
# of the terminal phase
# Half-life is multiplied by 3 and rounded to the closest 12 time units
  lambz <- -pred.lambdaz(obs[, 2], obs[, 1])[2]
  hl <- log(2)/lambz
  ceiling(hl/4)*12
}
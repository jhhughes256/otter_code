# Sum of exponential functions for sourcing
# -----------------------------------------------------------------------------
# The functions in order are:
#   pred.sumexp - gives sum of exponentials for given set of parameters
#   mle.sumexp - maximum likelihood estimation function to be minimised
#   optim.sumexp - determine parameters for differing numbers of exponentials
#   err.interv - trapezoidal error function to be minimised
#   optim.interv - determine intervals that give the smallest error
# -----------------------------------------------------------------------------
# Sum of exponentials predicted concentrations function
  # pred.sumexp <- function(par, x, d = 0) {
  # # There are currently issues with order.sumexp messing with optimisation if it is within pred.sumexp
  # # Provides predictions according to model parameters
  # # par = sum of exponential parameters
  # # x = independent variable (time)
  # # d = order of derivative (uses dth derivative of model)
  # # Define objects
  #   l <- length(par)  # number of parameters
  #   a <- l %% 2 == 1  # absorption status (odd parameter length == absorption)
  #   n <- ceiling(l/2)  # number of exponentials
  #   p <- order.sumexp(par, n, a)  # order parameters (allows for flip-flop)
  # # Sum of exponentials
  #   for (i in 1:n) {  # for each exponential
  #     if (i == 1) yhat <- p[i]^d*exp(p[i]*x + p[n+i])  # first exponential defines yhat
  #     else if (i != n | !a) yhat <- yhat + p[i]^d*exp(p[i]*x + p[n+i])  # following exponentials add to yhat
  #     else if (a) yhat <- yhat - p[i]^d*exp(p[i]*x)*sum(exp(p[(n+1):(2*n-1)]))  # for absorption curve apply final term
  #   }
  #   return(yhat)  # predicted dependent variable (drug concentration)
  # }

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

# Oldest function
# pred.sumexp <- function(x, t, d = 0) {
#   l <- length(x)
#   a <- ifelse(l %% 2 == 0, 0, 1)
#   n <- ceiling(l/2)
#   for (i in 1:n) {
#     if (i == 1) y <- x[i]^d*exp(x[i]*t + x[n+i])
#     else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
#     else if (a == 1) y <- y - x[i]^d*exp(x[i]*t)*sum(exp(x[(n+1):(2*n-1)]))
#   }
#   return(y)
# }

# Less Old function
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

# Maximum likelihood estimation function for parameter optimisation
  mle.sumexp <- function(par, x, y, sigma, ga = F) {
  # Determine log likelihood of given model parameters
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


# Fit sum of exponentials to curve for different numbers of exponentials
  # without hessian
  optim.sumexp <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    lmres <- unname(lm(log(y[lm.sub]) ~ x[lm.sub])$coefficients)
    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        optres <- list(
          par = c(lmres[2], lmres[1]),
          value = mle.sumexp(unname(c(lmres[2], lmres[1])), x, y, 0.01),
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        gares <- ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
          popSize = 250,
          monitor = F
        )
        optres <- optim(
          gares@solution[1, ],
          mle.sumexp,
          method = "BFGS",
          x = x, y = y, sigma = 0.01
        )
      }
      slope.par <- optres$par[1:(i+oral)]
      slope.ord <- order(slope.par, decreasing = T)
      par.ord <- unname(c(slope.par[slope.ord], optres$par[(i+oral+1):length(optres$par)]))
      opt.par[[i]] <- par.ord
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
      opt.mes[[i]] <- ifelse(is.null(optres$message),
        "NULL", optres$message)
    }
    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con, message = opt.mes)
    res
  }

  # with hessian matrix
# Fit sum of exponential parameters to data
  optim.sumexp.hes <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
  # Set up objects
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    opt.hes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    # browser()
    repeat {
      lmres <- unname(lm(log(y[lm.sub]) ~ x[lm.sub])$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -0.00001)
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- lm.sub[-1]
    }
    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        optres <- list(
          par = c(lmres[2], lmres[1]),
          value = mle.sumexp(unname(c(lmres[2], lmres[1])), x, y, 0.01),
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        gares <- try(ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
          popSize = 250,
          monitor = F
        ))
        if (class(gares) == "try-error") browser()
        optres <- try(optim(
          gares@solution[1, ],
          mle.sumexp,
          method = "BFGS", hessian = T,
          x = x, y = y, sigma = 0.01
        ))
        if (class(optres) == "try-error") {
          optres <- list(
            par = gares@solution[1,],
            value = mle.sumexp(gares@solution[1,], x, y, 0.01),
            counts = c("function" = 501, gradient = NA),
            convergence = 99,
            message = "zero gradient",
            hessian = matrix(NA,
              ncol = length(gares@solution[1,]),
              nrow = length(gares@solution[1,])
            )
          )
        }
      }
      slope.par <- -abs(head(optres$par, i + oral))
      int.par <- tail(optres$par, -(i + oral))
      slope.ord <- order(slope.par, decreasing = T)
      int.ord <- slope.ord  # intercept order (match slopes)
      if (oral) int.ord <- order(slope.par[-slope.ord[i]], decreasing = T)
      par.ord <- unname(c(slope.par[slope.ord], int.par[int.ord]))
      opt.par[[i]] <- par.ord
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
      opt.mes[[i]] <- ifelse(is.null(optres$message),
        "NULL", optres$message)
      opt.hes[[i]] <- optres$hessian
    }
    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con, message = opt.mes, hessian = opt.hes)
    res
  }

  optim.sumexp.se <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    opt.hes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    lmres <- unname(lm(log(y[lm.sub]) ~ x[lm.sub])$coefficients)
    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        optres <- list(
          par = c(lmres[2], lmres[1]),
          value = mle.sumexp(unname(c(lmres[2], lmres[1])), x, y, 0.01),
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        if (is.na(lmres[2])) {
          lmres <- c(max(y), -0.00001)
        }
        gares <- ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
          popSize = 250,
          monitor = F
        )
        optres <- optim(
          gares@solution[1, ],
          mle.sumexp,
          method = "BFGS", hessian = T,
          x = x, y = y, sigma = 0.01
        )
        repeat {
          vc_mat <- suppressWarnings(try(solve(optres$hessian)))
          if(class(vc_mat) != "try-error") {
            se <- suppressWarnings(sqrt(diag(vc_mat)))
            if (!any(is.nan(se))) {
              se_percent <- se/optres$par*100
              if (max(se_percent) <= 100) {
                break
              }
            }
          }
          gares <- ga("real-valued",
            mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
            min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
            max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
            selection = gareal_lrSelection,
            crossover = gareal_spCrossover,
            mutation = gareal_raMutation,
            maxiter = 50,
            popSize = 250,
            monitor = F
          )
          optres <- optim(
            gares@solution[1, ],
            mle.sumexp,
            method = "BFGS", hessian = T,
            x = x, y = y, sigma = 0.01
          )
        }
      }
      slope.par <- optres$par[1:(i+oral)]
      slope.ord <- order(slope.par, decreasing = T)
      par.ord <- unname(c(slope.par[slope.ord], optres$par[(i+oral+1):length(optres$par)]))
      opt.par[[i]] <- par.ord
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
      opt.mes[[i]] <- ifelse(is.null(optres$message),
        "NULL", optres$message)
      opt.hes[[i]] <- optres$hessian
    }
    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con, message = opt.mes, hessian = opt.hes)
    res
  }

  mle.sumexp.sig <- function(par, x, y, errmod, ga = F) {
  # Determine log likelihood of given model parameters
  # par = sum of exponential parameters
  # x = independent variable (time)
  # y = observed dependent variable (drug concentration)
  # errmod = error model to be used c("add", "prop", "both")
  # ga = genetic algorithm status
    nerr <- ifelse(errmod == "both", 2, 1) # set number of error parameters
    fit.par <- head(par, -nerr)  # define sum of exponentail parameters
    err.par <- tail(par, nerr)  # define error paramters
    z <- ifelse(ga, 2, -2)  # adjust ofv for minimisation or maximisation
    yhat <- pred.sumexp(fit.par, x)  # sum of exponential model prediction
  # Define standard deviation of normal distribution
    if (errmod == "add") sd <- err.par
    else if (errmod == "prop") sd <- abs(yhat)*err.par
    else if (errmod == "both") {
      add <- err.par[1]
      prop <- abs(yhat)*err.par[2]
      sd <- sqrt(add^2 + prop^2)
    }
    else stop("Please enter valid error model; \"add\", \"prop\" or \"both\"")
  # Determine log likelihood
    loglik <- suppressWarnings(dnorm(y, yhat, sd, log = T))
    return(z*sum(loglik))
  }

  optim.sumexp.sig <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
      error = NULL, hessian = NULL, message = NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    repeat {
      lm.mod <- lm(log(y[lm.sub]) ~ x[lm.sub])
      lmres <- unname(lm.mod$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -pred.lambdaz(y, x))
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- tail(lm.sub, -1)
    }
  # Estimate candidate model parameters
    lm.sd <- sd(residuals(lm.mod))
    lm.add <- matrix(c(0, max(y)*lm.sd), nrow = 2)
    lm.prop <- matrix(c(lm.sd/50, lm.sd*50), nrow = 2)
    cand.mod <- expand.grid(1:nexp, c("add", "prop", "both"))  # candidate models
    for (i in 1:nrow(cand.mod)) {
      mod <- cand.mod[i, ]
      if (mod[[2]] == "both") lm.err <- cbind(lm.add, lm.prop)
      else lm.err <- get(paste0("lm.", mod[[2]]))
      gares <- try(ga("real-valued",
        mle.sumexp.sig, x = x, y = y, ga = T, errmod = mod[[2]],
        min = c(rep(lmres[2]*50, mod[[1]] + oral), rep(lmres[1]-2, mod[[1]]), lm.err[1,]),
        max = c(rep(lmres[2]/50, mod[[1]] + oral), rep(lmres[1]+2, mod[[1]]), lm.err[2,]),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      optres <- try(optim(
        gares@solution[1, ],
        mle.sumexp.sig,
        method = "BFGS", hessian = T,
        x = x, y = y, errmod = mod[[2]]
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = mle.sumexp.sig(gares@solution[1,], x, y, errmod = mod[[2]]),
          counts = c("function" = 501, gradient = NA),
          convergence = 99,
          message = "zero gradient",
          hessian = matrix(NA,
            ncol = length(gares@solution[1,]),
            nrow = length(gares@solution[1,])
          )
        )
      }
    # Create output
      fit.par <- head(optres$par, -ncol(lm.err))
      err.par <- tail(optres$par, ncol(lm.err))
      par.ord <- order.sumexp(fit.par, mod[[1]], oral)
      res$par[[i]] <- optres$par
      res$sumexp[[i]] <- par.ord
      res$value[[i]] <- c(ofv = optres$value, optres$counts)
      res$error[[i]] <- c(signif(err.par, 5), type = paste(mod[[2]]))
      res$hessian[[i]] <- optres$hessian
      res$message[[i]] <- c(convergence = optres$convergence,
        message = ifelse(is.null(optres$message), "NULL", optres$message)
      )
    }
    res
  }

  mle.sumexp.err <- function(par, x, y, errmod, ga = F) {
  # Determine log likelihood of given model parameters
  # par = sum of exponential parameters
  # x = independent variable (time)
  # y = observed dependent variable (drug concentration)
  # errmod = error model to be used c("add", "prop", "both")
  # ga = genetic algorithm status
    nerr <- ifelse(errmod == "both", 2, 1) # set number of error parameters
    fit.par <- head(par, -nerr)  # define sum of exponentail parameters
    err.par <- tail(par, nerr)  # define error paramters
    z <- ifelse(ga, 2, -2)  # adjust ofv for minimisation or maximisation
    yhat <- pred.sumexp(fit.par, x)  # sum of exponential model prediction
  # Define standard deviation of normal distribution
    if (errmod == "add") sd <- err.par
    else if (errmod == "prop") sd <- abs(yhat)*err.par
    else if (errmod == "both") {
      add <- err.par[1]
      prop <- abs(yhat)*err.par[2]
      sd <- sqrt(add^2 + prop^2)
    }
    else stop("Please enter valid error model; \"add\", \"prop\" or \"both\"")
  # Determine log likelihood
    loglik <- suppressWarnings(dnorm(y, yhat, sd, log = T))
    return(z*sum(loglik))
  }

  optim.sumexp.new <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
    error = NULL, hessian = NULL, message = NULL)
    # opt.par <- list(NULL)
    # opt.val <- list(NULL)
    # opt.gra <- list(NULL)
    # opt.con <- list(NULL)
    # opt.mes <- list(NULL)
    # opt.hes <- list(NULL)
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    lmres <- -pred.lambdaz(y, x)[1]
    lm.sub <- which(y == max(y))[1]:length(y)
    repeat {
      lm.mod <- lm(log(y[lm.sub]) ~ x[lm.sub])
      lmres <- unname(lm.mod$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -pred.lambdaz(y, x))
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- tail(lm.sub, -1)
    }
  # Estimate candidate model parameters
    for (i in 1:nexp) {
      gares <- try(ga("real-valued",
        mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
        min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
        max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      optres <- try(optim(
        gares@solution[1, ],
        mle.sumexp,
        method = "BFGS", hessian = T,
        x = x, y = y, sigma = 0.01
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = mle.sumexp(gares@solution[1,], x, y, 0.01),
          counts = c("function" = 501, gradient = NA),
          convergence = 99,
          message = "zero gradient",
          hessian = matrix(NA,
            ncol = length(gares@solution[1,]),
            nrow = length(gares@solution[1,])
          )
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

# Chi-squared difference test
# Takes a list of optim results and gives the best optim result
  chisq.sumexp <- function(opt) {
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- opt$value[[i]] - opt$value[[j]]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.01) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }

  chisq.sumexp.aic <- function(opt) {
    x <- unlist(opt$value)
    k <- unlist(lapply(opt$par, length))
    aic <- x + 2*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

  chisq.sumexp.bic <- function(opt, nobs) {
    x <- unlist(opt$value)
    k <- unlist(lapply(opt$par, length))
    bic <- x + log(nobs)*k
    return(sapply(opt, function(x) x[which(bic == min(bic))]))
  }

  best.sumexp.lrt <- function(opt) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- ofv[i] - ofv[j]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.01) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }

  best.sumexp.aic <- function(opt) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    k <- unlist(lapply(opt$par, length))
    aic <- ofv + 2*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

  best.sumexp.bic <- function(opt, nobs) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    k <- unlist(lapply(opt$par, length))
    bic <- ofv + log(nobs)*k
    return(sapply(opt, function(x) x[which(bic == min(bic))]))
  }

# Trapezoidal error function for interval optimisation
  # standard using times
  err.interv <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
    times <- c(tfirst, par, tlast, tmax)
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    sum(err^2)
  }

  # ga using times
  err.interv.ga <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F, ga = F) {
    z <- ifelse(ga, -1, 1)
    times <- unique(c(tfirst, par, tlast, tmax))
    times <- times[order(times)]
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    return(z*sum(err^2))
  }

  # standard using dt
  err.interv.dt <- function(par, exp.par, tfirst, tlast, a = F) {
    times <- c(tfirst, cumsum(par), tlast)
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    return(sum(err^2))
  }

  err.interv.dtmax <- function(par, exp.par, tfirst, tlast, a = F, tmax = NULL) {
    times <- c(tfirst, cumsum(par), tmax, tlast)
    times <- sort(times)
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    sumerr <- sqrt(mean(err^2))
    return(sumerr)
  }

# Interval optimising function
  # standard using times
  optim.interv <- function(times, par, tmax = NULL) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
    if (!is.null(tmax)) init.par <- init.par[-length(init.par)]
    xmin <- min(x)
    xmax <- max(x)
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    res <- optim(
      init.par,
      err.interv,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = par, tfirst = xmin + 0.01, tlast = xmax - 0.01, tmax = tmax,
      a = absorp
    )
    return(res)
  }

  # using dt instead of times
  optim.interv.dt <- function(par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    npar <- length(times) - 2
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    init.par <- cumsum(rep(tlast/48, npar))
    res <- optim(
      init.par,
      err.interv.dt,
      method = "L-BFGS-B", hessian = T,
      lower = tlast/48, upper = tlast - npar*tlast/48,
      exp.par = par, tfirst = tfirst, tlast = tlast, a = absorp
    )
    res$times <- cumsum(res$par)
    return(res)
  }

  optim.interv.dtmax <- function(par, times, tmax = FALSE) {
    tfirst <- min(times)
    tlast <- max(times)
    tbound <- tlast - tfirst
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    npar <- length(times) - (2 + tmax*absorp)
    init.par <- cumsum(rep(tbound/48, npar))
    if (tmax & absorp) tmax.val <- tmax.sumexp(par, tlast, tlast/720)  # maximum length of 721
    else tmax.val <- NULL
    res <- try(optim(
      init.par,
      err.interv.dtmax,
      method = "L-BFGS-B", hessian = T,
      lower = tbound/48, upper = tbound/2, exp.par = par,
      tfirst = tfirst, tlast = tlast, a = absorp, tmax = tmax.val
    ))
    if (class(res) == "try-error") {
      res <- try(optim(
        init.par,
        err.interv.dtmax,
        method = "L-BFGS-B", hessian = T,
        lower = tbound/48, upper = tbound/(npar/2), exp.par = par,
        tfirst = tfirst, tlast = tlast, a = absorp, tmax = tmax.val
      ))
    }
    res$times <- sort(c(cumsum(c(tfirst, res$par)), tmax.val, tlast))
    return(res)
  }

  # using ga for initial parameters
  optim.interv.ga100 <- function(par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    if (is.null(tmax)) {
      is.tmax <- 2
    } else if (tmax == 0) {
      is.tmax <- 2
    } else {
      is.tmax <- 3
    }
    npar <- length(times) - is.tmax
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    flag <- 0
    repeat {
      gares <- ga("real-valued",
        err.interv.ga, exp.par = par, a = absorp,
        tfirst = tfirst, tlast = tlast, tmax = tmax, ga = T,
        min = rep(tfirst + 0.01, npar), max = rep(tlast - 0.01, npar),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      )
      res <- optim(
        gares@solution[order(gares@solution)],
        err.interv.ga,
        method = "BFGS", hessian = T,
        exp.par = par, tfirst = tfirst, tlast = tlast, tmax = tmax,
        a = absorp
      )
      vc_mat <- try(solve(res$hessian))
      if(class(vc_mat) != "try-error") {
        se <- sqrt(diag(vc_mat))
        if (!any(is.nan(se))) {
          se_percent <- se/res$par*100
          if (max(se_percent) <= 100) {
            res$flag <- flag
            res$se <- se_percent
            break
          }
        }
      }
      flag <- flag + 1
      if (flag == 10) {
        res$flag <- flag
        res$se <- NA
        break
      }
    }
    return(res)
  }

  optim.interv.ga50 <- function(par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    is.tmax <- ifelse(is.null(tmax), 2, 3)
    npar <- length(times) - is.tmax
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    flag <- 0
    repeat {
      gares <- ga("real-valued",
        err.interv.ga, exp.par = par, a = absorp,
        tfirst = tfirst, tlast = tlast, tmax = tmax, ga = T,
        min = rep(tfirst + 0.01, npar), max = rep(tlast - 0.01, npar),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      )
      res <- optim(
        gares@solution[order(gares@solution)],
        err.interv.ga,
        method = "BFGS", hessian = T,
        exp.par = par, tfirst = tfirst, tlast = tlast, tmax = tmax,
        a = absorp
      )
      vc_mat <- try(solve(res$hessian))
      if(class(vc_mat) != "try-error") {
        se <- sqrt(diag(vc_mat))
        if (!any(is.nan(se))) {
          se_percent <- se/res$par*100
          if (max(se_percent) <= 50 | flag == 10) {
            res$flag <- flag
            res$se <- se_percent
            break
          }
        }
      }
      flag <- flag + 1
    }
    return(res)
  }

  pred.tlast <- function(par, tlast) {
    i <- round(tlast/12, 0)
    perc.term <- 1
    timeloop <- seq(0, i*12, by = i*12/120)
    predloop <- pred.sumexp(par, timeloop)
    tmax <- timeloop[which(predloop == max(predloop))]
    while(perc.term > 0.2 & i < 730) {
      if (exists("init")) {
        repeat {
          i <- i + 1
          timeloop <- seq(0, i*12, by = i*12/120)
          predloop <- pred.sumexp(par, timeloop)
          clast <- tail(predloop, 1)
          if (clast < cterm | i == 730) break
        }
      }
      clast <- tail(predloop, 1)
      auclast <- auc.interv.sumexp(timeloop, par, log = T)
      lambz <- max(head(par, ceiling(length(par)/2)))
      aucinf <- clast/-lambz
      perc.term <- aucinf/(auclast+aucinf)
      cterm <- clast*(0.18/perc.term)
      init <- 1
    }
    return(c(i*12, 1-perc.term))
  }

  pred.tlast.lam <- function(par) {
    nexp <- ceiling(length(par)/2)
    lambz <- abs(max(par[1:nexp]))
    hl <- log(2)/lambz
    ceiling(hl/4)*12
  }

  obs.tlast.lam <- function(obs) {
    lambz <- pred.lambdaz(obs[, 2], obs[, 1])
    hl <- log(2)/lambz
    ceiling(hl/4)*12
  }

  pred.lambdaz <- function(dv, t) {
    if (t[1] == 0) dv[1] <- 0
    mdv <- which(dv == 0)
    i <- 3
    bestR2 <- -1
    bestk <- 0
    if (length(dv[-mdv]) >= i) {
      repeat {
        mod <- suppressWarnings(lm(log(tail(dv[-mdv], i)) ~ tail(unique(t)[-mdv], i)))
        k <- -1*mod$coefficients["tail(unique(t)[-mdv], i)"]
        R2 <- suppressWarnings(as.numeric(summary(mod)["adj.r.squared"]))
        if (is.na(k)) browser()
        if (is.nan(R2)) R2 <- suppressWarnings(as.numeric(summary(mod)["r.squared"]))
        if (k > 0) {
          if (R2 > bestR2) {
            if (i > 2) bestR2 <- R2
            bestk <- k
          } else {
            break
          }
        }
        if (i == 5 & bestk == 0) {  #
          mdv <- c(mdv, which(dv == max(tail(dv[-mdv], 3))))
          i <- 1
        }
        if (i == length(dv[-mdv])) {  # last ditch effort (intended for simulation study only)
          if (bestk > 0) break
          else {
            mod <- suppressWarnings(lm(log(tail(dv, 2)) ~ tail(unique(t), 2)))
            bestk <- -1*mod$coefficients["tail(unique(t), 2)"]
            if (bestk > 0) break
            else {
              bestk <- log(2)/56
              break
            }
          }
        }
        i <- i + 1
      }
    }
    bestk
  }

# -----------------------------------------------------------------------------
# Determine tmax given a set of sumexp parameters
  tmax.sumexp <- function(par, tlast = 24, res = 0.1) {
    times <- seq(0, tlast, by = res)
    yhat <- pred.sumexp(par, times)
    return(times[which(yhat == max(yhat))])
  }

# Determine auc given a set of intervals
  auc.interv <- function(times, par, fn, log = F) {
    C <- do.call(fn, list(x = times, p = par))
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }

  auc.interv.sumexp <- function(times, par, log = F) {
    C <- pred.sumexp(par, times)
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }

  auc.interv.lam <- function(par, t) {
    dv <- pred.sumexp(par, unique(t))
    lambz <- pred.lambdaz(dv, unique(t))
    tail(dv, 1)/lambz
  }

# Without for loop
#  auc.interv <- function(times, par, fn, log = F) {
#    C <- do.call(fn, list(x = times, p = par))
#    h <- diff(times)
#    EC <- sum(C[-1], C[-length(C)])
#    dC <- diff(-C)
#    if (!log) auc <- EC*h/2
#    else auc[i-1] <- dC*h/log(C[i-1]/C[i])
#    return(sum(auc))
#  }

# Plot random data
  plot.rdata <- function(data, t, n, interv) {
    plotdata <- data.frame(
      id = rep(1:n, each = length(t)),
      time = rep(t, times = n),
      dv = as.vector(data)
    )
    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
    plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
    plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)")
    plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = "free")
    return(plotobj)
  }

# Random sum of exponential data for sourcing
# -----------------------------------------------------------------------------
# The datasets in order are:
#   d1b - one compartment kinetics given iv
#   d2b - two compartment kinetics given iv
#   d3b - three compartment kinetics given iv
#   d1a - one compartment kinetics given orally
#   d2a - two compartment kinetics given orally
#   d3a - three compartment kinetics given orally
# -----------------------------------------------------------------------------
# Uncomment below code to view random data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# library(ggplot2)
# theme_bw2 <- theme_set(theme_bw(base_size = 14))
# theme_update(plot.title = element_text(hjust = 0.5))
# plot.rdata <- function(data, t, n, interv, log = F) {
#   plotdata <- data.frame(
#     id = rep(1:n, each = length(t)),
#     time = rep(t, times = n),
#     dv = as.vector(data)
#   )
#   xlim <- c(t[1], t[length(t)])
#   plotobj <- NULL
#   plotobj <- ggplot(data = plotdata)
#   plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
#   plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
#   plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
#   if (!log) plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
#   else plotobj <- plotobj + scale_y_log10("Concentration (mg/mL)\n")
#   plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
#   plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = "free")
#   return(plotobj)
# }
# -----------------------------------------------------------------------------
  if (!exists("niter")) niter <- 16
  time.samp <- seq(0, 24, by = 0.1)
# One Compartment Kinetics
  pred.d1b <- function(x, p) {
    exp(p[1]*x + p[2])
  }
  d1b.p <- matrix(nrow = 2, ncol = niter)
  d1b.p[1, ] <- runif(niter, -0.8, -0.08)
  d1b.p[2, ] <- runif(niter, 1, 6)
  d1b <- apply(d1b.p, 2, function(p, x) pred.d1b(x, p), x = time.samp)
  d1b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)
  #plot.rdata(d1b, time.samp, niter, -1, log = F)

# Two Compartment Kinetics
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  d2b.p <- matrix(nrow = 4, ncol = niter)
  d2b.p[1, ] <- runif(niter, -1, -0.1)
  d2b.p[2, ] <- runif(niter, d2b.p[1, ]*0.8, d2b.p[1, ]*0.05)
  d2b.p[3, ] <- runif(niter, 1, 6)
  d2b.p[4, ] <- runif(niter, d2b.p[3, ] - 1, d2b.p[3, ] - 0.2)
  d2b <- apply(d2b.p, 2, function(p, x) pred.d2b(x, p), x = time.samp)
  d2b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)
  #plot.rdata(d2b, time.samp, niter, -1, log = F)

# Three Compartment Kinetics
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  d3b.p <- matrix(nrow = 6, ncol = niter)
  d3b.p[1, ] <- runif(niter, -1, -0.1)
  d3b.p[2, ] <- runif(niter, d3b.p[1, ]*0.8, d3b.p[1, ]*0.05)
  d3b.p[3, ] <- runif(niter, d3b.p[2, ]*0.8, d3b.p[2, ]*0.05)
  d3b.p[4, ] <- runif(niter, 1, 6)
  d3b.p[5, ] <- runif(niter, d3b.p[4, ] - 1, d3b.p[4, ] - 0.2)
  d3b.p[6, ] <- runif(niter, d3b.p[5, ] - 1, d3b.p[5, ] - 0.2)
  d3b <- apply(d3b.p, 2, function(p, x) pred.d3b(x, p), x = time.samp)
  d3b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)
  #plot.rdata(d3b, time.samp, niter, -1, log = F)

# One Compartment Kinetics w/ Absorption
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.p <- matrix(nrow = 3, ncol = niter)
  d1a.sim <- 1:niter
  repeat {
    nsim <- length(d1a.sim)
    d1a.p[1, d1a.sim] <- runif(nsim, -1, -0.1)
    d1a.p[2, d1a.sim] <- runif(nsim, d1a.p[1, d1a.sim]*2, d1a.p[1, d1a.sim]*1.1)
    d1a.p[3, d1a.sim] <- runif(nsim, 1, 6)
    d1a <- apply(d1a.p, 2, function(p, x) pred.d1a(x, p), x = time.samp)
    d1a.tmax <- apply(d1a, 2, function(x) time.samp[which(x == max(x))])
    d1a.sim <- which(d1a.tmax > 12)
    if (length(d1a.sim) == 0) break
    else print(paste("d1a repeat", length(d1a.sim)))
  }
  d1a.t <- c(0, 0.5, 1, 3, 5, 7, 10, 16, 24)
  #plot.rdata(d1a, time.samp, niter, -1, log = F)

# Two Compartment Kinetics w/ Absorption
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  d2a.p <- matrix(nrow = 5, ncol = niter)
  d2a.sim <- 1:niter
  repeat {
    nsim <- length(d2a.sim)
    d2a.p[1, d2a.sim] <- runif(nsim, -1, -0.1)
    d2a.p[2, d2a.sim] <- runif(nsim, d2a.p[1, d2a.sim]*0.8, d2a.p[1, d2a.sim]*0.05)
    d2a.p[3, d2a.sim] <- runif(nsim, d2a.p[1, d2a.sim]*2, d2a.p[1, d2a.sim]*1.1)
    d2a.p[4, d2a.sim] <- runif(nsim, 1, 6)
    d2a.p[5, d2a.sim] <- runif(nsim, d2a.p[4, d2a.sim] - 1, d2a.p[4, d2a.sim] - 0.2)
    d2a <- apply(d2a.p, 2, function(p, x) pred.d2a(x, p), x = time.samp)
    d2a.tmax <- apply(d2a, 2, function(x) time.samp[which(x == max(x))])
    d2a.sim <- which(d2a.tmax > 12)
    if (length(d2a.sim) == 0) break
    else print(paste("d2a repeat", length(d2a.sim)))
  }
  d2a.t <- c(0, 0.5, 1, 3, 5, 8, 12, 16, 24)
  #plot.rdata(d2a, time.samp, niter, -1, log = F)

# Three Compartment Kinetics w/ Absorption
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  d3a.p <- matrix(nrow = 7, ncol = niter)
  d3a.sim <- 1:niter
  repeat {
    nsim <- length(d3a.sim)
    d3a.p[1, d3a.sim] <- runif(nsim, -1, -0.1)
    d3a.p[2, d3a.sim] <- runif(nsim, d3a.p[1, d3a.sim]*0.8, d3a.p[1, d3a.sim]*0.05)
    d3a.p[3, d3a.sim] <- runif(nsim, d3a.p[2, d3a.sim]*0.8, d3a.p[2, d3a.sim]*0.05)
    d3a.p[4, d3a.sim] <- runif(nsim, d3a.p[1, d3a.sim]*2, d3a.p[1, d3a.sim]*1.1)
    d3a.p[5, d3a.sim] <- runif(nsim, 1, 6)
    d3a.p[6, d3a.sim] <- runif(nsim, d3a.p[5, d3a.sim] - 1, d3a.p[5, d3a.sim] - 0.2)
    d3a.p[7, d3a.sim] <- runif(nsim, d3a.p[6, d3a.sim] - 1, d3a.p[6, d3a.sim] - 0.2)
    d3a <- apply(d3a.p, 2, function(p, x) pred.d3a(x, p), x = time.samp)
    d3a.tmax <- apply(d3a, 2, function(x) time.samp[which(x == max(x))])
    d3a.sim <- which(d3a.tmax > 12)
    if (length(d3a.sim) == 0) break
    else print(paste("d3a repeat", length(d3a.sim)))
  }
  d3a.t <- c(0, 0.5, 1.5, 3, 5, 8, 12, 16, 24)
  #plot.rdata(d3a, time.samp, niter, -1, log = F)

# Sum of exponential data for sourcing
# -----------------------------------------------------------------------------
# The datasets in order are:
#   d1b - one compartment kinetics given iv
#   d2b - two compartment kinetics given iv
#   d3b - three compartment kinetics given iv
#   d1a - one compartment kinetics given orally
#   d2a - two compartment kinetics given orally
#   d3a - three compartment kinetics given orally
# -----------------------------------------------------------------------------
  time.samp <- seq(0, 24, by = 0.05)

# Two Compartment Kinetics
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  d2b.p <- c(-0.06, -0.5, 5, 6)
  d2b <- pred.d2b(time.samp, d2b.p)
  #with(d2b, plot(time, log(conc)))
  d2b.t <- c(0, 0.5, 1, 3, 5, 8, 12, 16, 24)

# Three Compartment Kinetics
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  d3b.p <- c(-0.03, -0.15, -0.5, 3.8, 4.7, 6)
  d3b <- pred.d3b(time.samp, d3b.p)
  #with(d3b, plot(time, log(conc)))
  d3b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)

# One Compartment Kinetics w/ Absorption
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.p <- c(-0.2, -0.4, 4)
  d1a <- pred.d1a(time.samp, d1a.p)
  #with(d1a plot(time, log(conc)))
  d1a.t <- c(0, 0.5, 1, 3, 5, 7, 10, 16, 24)

# Two Compartment Kinetics w/ Absorption
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  d2a.p <- c(-0.06, -0.4, -0.6, log(exp(4)*0.15), log(exp(4)*0.85))
  d2a <- pred.d2a(time.samp, d2a.p)
  # with(d2a, plot(time, log(conc)))
  d2a.t <- c(0, 0.5, 1, 3, 5, 8, 12, 16, 24)

# Three Compartment Kinetics w/ Absorption
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  d3a.p <- c(-0.05, -0.25, -0.6, -1, log(exp(4)*0.10), log(exp(4)*0.35), log(exp(4)*0.55))
  d3a <- pred.d3a(time.samp, d3a.p)
  #with(d3a, plot(time[-1], log(conc)[-1]))
  d3a.t <- c(0, 0.5, 1.5, 3, 5, 8, 12, 16, 24)

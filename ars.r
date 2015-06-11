# Derivative-free adaptive rejection sampling.
# Author: Urban Soban

ars <- function (f, domain, init.points=c(-4, 4), n.samples=100) {
  # Adaptive rejection sampling.
  # 
  # Args:
  #   f: log-concave density function to sample. Should accept only one parameter, at which it is evaluated.
  #   y: domain of support (e.g. c(-Inf, Inf))
  #   init.points: min and max theta to sample
  #   n.samples: number of samples to take
  #
  # Returns:
  #   Vector of samples taken from distribution f
  h <- function (x) { log(f(x)) } # log of f.

  # Choose starting (minimal 3) points on abscisa.
  start.points <- function () {
    if (is.infinite(domain[1]) && is.infinite(domain[2])) {
      pts <- c(init.points[1], runif(1, init.points[1], init.points[2]), init.points[2])
    } 
    else if (is.finite(domain[1])) {
      pts <- c(domain[1] + 0.0001, runif(1, domain[1] + 0.0001, init.points[2]), init.points[2])
    } 
    else {
      pts <- init.points
    }
    pts
  }
  
  is.logconcave <- function () {
    if (h.chords$k[1] > 0 && h.chords$k[-nrow(h.chords)] < 0) return(TRUE)
    else return(FALSE)
  }
  
  # Makes chord lines passing throught neighboring points lying on h=log(f)
  mk.chords <- function () {
    np <- nrow(h.points)
    ks <- (h.points$y[-np] - h.points$y[-1])/(h.points$x[-np] - h.points$x[-1])
    ns <- h.points$y[-np] - ks*h.points$x[-np]
    data.frame(k=ks, n=ns, xstart=h.points$x[-np], xend=h.points$x[-1])
  }
  
  # Makes upper hull.
  mk.uphull <- function () {
    n.pts <- nrow(h.points)
    cuts <- c()
    ks <- c()
    ns <- c()
    for (i in seq(1, nrow(h.points))) {
      cuts <- c(cuts, h.points$x[i])
      if (i == 1) {
        # Edge case 1. @TODO
        ks <- c(ks, h.chords$k[i+1])
        ns <- c(ns, h.chords$n[i+1])
      }
      else if (i == n.pts || i == n.pts - 1) {
        # Edge case 2. @TODO
        ks <- c(ks, h.chords$k[i-1])
        ns <- c(ns, h.chords$n[i-1])
      }
      else {
        # Compute intersection.
        isect.x <- (h.chords$n[i+1] - h.chords$n[i-1])/(h.chords$k[i-1] - h.chords$k[i+1])
        ks <- c(ks, h.chords$k[i-1], h.chords$k[i+1])
        ns <- c(ns, h.chords$n[i-1], h.chords$n[i+1])
        cuts <- c(cuts, isect.x)
      }
    }
    data.frame(cut=cuts, k=ks, n=ns)
  }
  
  # Computes proportions.
  mk.props <- function () {
    ks <- head(up.hull$k, -1)
    ns <- head(up.hull$n, -1)
    theta_curr <- head(up.hull$cut, -1)
    theta_nxt  <- tail(up.hull$cut, -1) 
    #browser()
    props <- exp(ns)*((exp(ks*theta_nxt) - exp(ks*theta_curr))/ks)
    props
  }
  
  evaluate.uphull <- function (theta) {
    params.idx <- findInterval(theta, up.hull$cut)
    up.hull$k[params.idx]*theta + up.hull$n[params.idx]
  }
  evaluate.lohull <- function (theta) {
    params.idx <- findInterval(theta, h.chords$xstart)
    h.chords$k[params.idx]*theta + h.chords$n[params.idx]
  }
  
  inv.cdf.sample <- function () {
    Gj <- sum(up.hull.props)
    probs <- up.hull.props/Gj
    s.i <- which.max(rmultinom(1, 1, probs))
    s.k <- up.hull$k[s.i]
    u <- runif(1)
    theta <- (1/s.k) * log(exp(s.k*up.hull$cut[s.i]) + u*(exp(s.k*up.hull$cut[s.i+1]) - exp(s.k*up.hull$cut[s.i])))
    theta
  }
  
  h.points <- start.points()
  h.points <- data.frame(x=h.points, y=sapply(h.points, h))
  h.chords <- mk.chords()
  up.hull  <- mk.uphull()
  up.hull.props <- mk.props()
  up.hull.norm <- sum(up.hull.props)

  n.sampled <- 0
  ars.samples = c()
  while (n.sampled < n.samples) {
    theta <- inv.cdf.sample()
    u <- runif(1)
    up.hull.theta <- evaluate.uphull(theta)
    lo.hull.theta <- evaluate.lohull(theta)
    if (u <= exp(lo.hull.theta - up.hull.theta)) {
      ars.samples <- c(ars.samples, theta)
      n.sampled <- n.sampled + 1
    }
    else {
      h.theta <- h(theta)
      if (u <= exp(h.theta-up.hull.theta)) {
        ars.samples <- c(ars.samples, theta)
        n.sampled <- n.sampled + 1
      }
      h.points <- rbind(h.points, data.frame(x=theta, y=h.theta))
      h.points <- h.points[with(h.points, order(x)),]
      h.chords <- mk.chords()
      up.hull <- mk.uphull()
      up.hull.props <- mk.props()
      up.hull.norm <- sum(up.hull.props)
    }
  }
  
  ars.samples
}
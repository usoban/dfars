# Rejection sampling.
rs <- function (f, prop.r, prop.d, n.samples=100) {
  n <- 0
  #rejected <- 0
  #tries <- 0
  #n.fails <- 0
  samples <- c()
  while (n < n.samples) {
    #tries <- tries + 1
    x <- prop.r(1)
    fx <- f(x)
    px <- prop.d(x)
    u <- runif(1)
    
    if (u < fx/px) {
      samples <- c(samples, x)
      n <- n + 1
      #rs.ntries <<- c(rs.ntries, n.fails)
      #n.fails <- 0
    }
    #else {
    #  rejected <- rejected + 1
    #  n.fails <- n.fails + 1
    #}
  }
  #rs.rejections <<- rejected
  #rs.tries <<- tries
  samples
}

plt.i <- 0
ars_plot <- function (f, x.min, x.max, pts, lo.hull, up.hull) {
  library(ggplot2)
  lgf <- function (x) {
    log(f(x))
  }
  x <- seq(x.min, x.max, by=0.01)
  y <- sapply(x, lgf)
  plt <- data.frame(x=x, y=y)
  #browser()
  pts.lo.hull <- data.frame(x=head(pts$x, -1), y=head(pts$y, -1), xend=tail(pts$x, -1), yend=tail(pts$y, -1), colour=rep('h', length(pts$y)-1))
  pts.hi.hull.x = c()
  pts.hi.hull.y = c()
  pts.hi.hull.endx = c()
  pts.hi.hull.endy = c()
  for (i in seq(1, nrow(up.hull))) {
    pts.hi.hull.x <- c(pts.hi.hull.x, up.hull$cut[i])
    pts.hi.hull.y <- c(pts.hi.hull.y, up.hull$k[i]*up.hull$cut[i] + up.hull$n[i])
    pts.hi.hull.endx <- c(pts.hi.hull.endx, up.hull$cut[i+1])
    pts.hi.hull.endy <- c(pts.hi.hull.endy, up.hull$k[i]*up.hull$cut[i+1] + up.hull$n[i])
  }
  pts.hi.hull <- data.frame(x=pts.hi.hull.x, y=pts.hi.hull.y, xend=pts.hi.hull.endx, yend=pts.hi.hull.endy)
  g1 <- ggplot(NULL, aes(x=x, y=y)) +
    geom_line(data=plt, mapping=aes(x=x, y=y, colour="h(theta)"), size=1.5) +
    geom_segment(data=pts.lo.hull, mapping=aes(x=x, y=y, xend=xend, yend=yend, colour='lower hull'), size=1) +
    geom_segment(data=pts.hi.hull, mapping=aes(x=x, y=y, xend=xend, yend=yend, colour='upper hull'), size=1) +
    geom_point(data=pts, mapping=aes(x=x, y=y)) + theme(legend.position='none')
  #guides(colour=guide_legend(title='', direction='horizontal', label.position='bottom')) +
  #theme(legend.position='bottom') + xlab(expression(theta)) + ylab('log f(theta)') + ggtitle('N(0, 0.2)')
  
  #   g2 <- ggplot(NULL, aes(x=x, y=y)) +
  #           geom_line(data=plt, mapping=aes(x=x, y=y), size=1.5, color='#CC5B75') +
  #           geom_segment(data=pts.lo.hull, mapping=aes(x=x, y=exp(y), xend=xend, yend=exp(yend)), size=1, color='#97DD6D') +
  #           geom_segment(data=pts.hi.hull, mapping=aes(x=x, y=y, xend=xend, yend=yend), size=1, color='#7D5CAB') +
  #           geom_point(data=pts, mapping=aes(x=x, y=y)) #+
  #scale_y_continuous(limits = c(-25, 25))
  
  ggsave(paste('/home/usoban/workspace/fri-mag/fri-math/bayes/ars/anim/plt-', plt.i, '.pdf', sep=''), plot=g1, width = 6, height=3)
  plt.i <<- plt.i + 1
}

# TESTS 
cmp.plt <- function (n, dist1, dist2, fdens, xmin=0, xmax=1) {
  samples1 <- dist1(n)
  xs <- seq(xmin+(1/n), xmax, by=(1/n))
  tmd.s <- proc.time()
  dens <- fdens(xs)
  tmd.e <- proc.time()
  print(tmd.e - tmd.s)
  dta <- data.frame(dist1=samples1, dist2=dist2, densx=xs, densy=dens)
  print(ggplot(dta) + geom_density(aes(x=dist1, colour='rbeta(2, 5)')) + geom_density(aes(x=dist2, colour='aes(beta(2, 5))')) + geom_line(aes(x=densx, y=densy, colour='beta(2, 5)')))
}
cmp.plt.nonprop <- function (n, dist1, fdens, xmin=0, xmax=1) {
  xs <- seq(xmin+(1/n), xmax, by=(1/n))
  dens <- fdens(xs)
  d1max <- max(density(dist1)$y)
  d2max <- max(dens)
  coef <- d1max/d2max
  dta <- data.frame(dist1=dist1, densx=xs, densy=coef*dens)
  g1 <- ggplot(dta) + 
    geom_density(aes(x=dist1, colour='ARS p(Beta|y, x, Alpha)'), size=1.3) + 
    geom_line(aes(x=densx, y=densy, colour='p(Beta|y, x, Alpha)'), size=1.2) +
    #xlab(expression(beta)) + ggtitle(paste('Sampling p(beta|x, y, alpha). Rejected: ', ars.rejections, '(', ars.rejections/10, '%)', split='')) + 
    guides(colour=guide_legend(title='')) + theme(legend.position='bottom')
  print(g1)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## BETA
# beta_dist <- function (alpha, beta) {
#   b <- beta(alpha, beta)
#   function (x) {
#     (x^(alpha-1)*(1-x)^(beta-1))/b
#   }
# }
# beta_rand <- function (alpha, beta) {
#   function (n) {
#     rbeta(n, shape1=alpha, shape2=beta)
#   }
# }

# tma.s <- proc.time()
# beta.samples <- ars(beta_dist(2, 5), domain=c(0, 1), init.points=c(0.0001, 0.9999), 10000)
# tma.e <- proc.time()
# print(tma.e - tma.s)
# cmp.plt(10000, beta_rand(2, 5), beta.samples, fdens=function(xs) { dbeta(xs, shape1=2, shape2=5) })

## GAMMA
# gamma_dist <- function (shape, scale) {
#   g <- gamma(shape)
#   function (x) {
#     (x^(shape-1) * exp(-x/scale))/(g*(scale^shape))
#   }
# }
# gamma_rand <- function (shape, scale) {
#   function (n) {
#     rgamma(n, shape=shape, scale=scale)
#   }
# }
# 
# gamma.samples <- ars(gamma_dist(3, 2), domain=c(0, Inf), init.points=c(0.0001, 20), n.samples=1000)
# gd <- gauss_dist(2, 3)
# gamma.samples.rs <- rs(gamma_dist(2, 2), gauss_rand(2, 3), function(x){ 3*gd(x) }, n.samples=1000)
# cmp.plt(1000, gamma_rand(2, 2), gamma.samples.rs, fdens=function(xs) { dgamma(xs, shape = 2, scale = 2) }, xmin=0, xmax=20)
# 
# print(paste("Rejects RS:", rs.rejections/rs.tries, "Rejects ARS:", ars.rejections/ars.tries))

#fails.dta <- data.frame(x=1:1000, yrs=rs.ntries, yars=ars.ntries)
#fail.g <- ggplot(fails.dta) + geom_line(aes(x=x, y=yrs, colour='Rejection sampling: 3*N(2, 3) (Rejected: 67.3%)'), size=1.3) + 
#                              geom_line(aes(x=x, y=yars, colour='Adaptive rejection sampling (Rejected: 1.7%)'), size=1.3) +
#                              xlab('Sample') + ylab('N draws') + theme(legend.position='bottom') +
#                              guides(color=guide_legend(title='')) + ggtitle('Number of draws per sample - G(2, 2)')
#print(fail.g)

## GAUSSIAN
# gauss_dist <- function (mu, sigma) {
#   function (x) {
#     exp(-(x-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi))
#   }
# }
# gauss_rand <- function (mu, sigma) {
#   function (n) {
#     rnorm(n, mu, sigma)
#   }
# }
# 
# gauss.samples <- ars(gauss_dist(0, 0.2), domain=c(-Inf, Inf), init.points=c(-2, 2), n.samples=1000)
# cmp.plt(1000, gauss_rand(0, 0.2), gauss.samples, fdens=function(xs) { dnorm(xs, 0, 0.2) }, xmin=-2, xmax=2)

# LOGIT
logit <- function(beta, 
                  y = c(0,1,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1), 
                  x=c(13.03,13.69,12.62,11.70,12.39,12.44,12.22,13.65,14.30,12.39,13.51,13.21,13.02,12.89,14.79,11.57,14.26,13.94,13.85,11.34,12.41,13.88,11.98,11.18,13.17,14.62,15.41,14.70,13.45,12.66,12.73,12.36,10.59,12.02,13.52,12.97,12.96,11.42,13.12,13.37,13.11,13.10,11.83),
                  alpha=0.01273263,
                  mu0 = 0, 
                  sigma20 = 1000)
{
  lik = dnorm(beta, mu0, sqrt(sigma20), log = T)
  for (i in 1:length(y))
  {
    p.i = 1 / (1 + exp(-alpha - beta*x[i]))
    lik = lik + dbinom(y[i], 1, p.i, log = T)
  }
  exp(lik)
}

logit.dens <- function () {
  function (xs) {
    logit(xs) 
  }
}

logit.samples <- ars(logit.dens(), domain=c(-Inf, Inf), init.points=c(-.1, .15), n.samples=10000)
cmp.plt.nonprop(10000, logit.samples, fdens=logit.dens(), xmin=-.1, xmax=.15)

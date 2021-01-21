library(mixdist)
library(gridExtra)

header <- c('bin.number', 'bin.diameter', 'perc.vol', 'v4', 'v5', 'v6')


read.coulter.data <- function(fpath) {
  data <- read.csv(fpath, header=FALSE, skip=55)
  names(data) <- header
  data <- data[complete.cases(data),]
  data$bin.max = c(data$bin.diameter[-1], max(data$bin.diameter) + 0.350)
  binwidths = data$bin.max - data$bin.diameter
  data$percvol.density = data$perc.vol/binwidths
  
  data
}

find.normal.fit <- function(data) {
  fit.normal <- nls(data$percvol.density ~ 
                      C1 * exp(-(data$bin.diameter-mean1)**2/(2 * sigma1**2)) + 
                      C2 * exp(-(data$bin.diameter-mean2)**2/(2 * sigma2**2)),
                    data=data,
                    start=list(C1=5, mean1=7, sigma1=2.2, C2=5, mean2=21, sigma2=5.2))
  
  fitted.normal.values <- coef(fit.normal)
  C1 <- fitted.normal.values['C1']
  C2 <- fitted.normal.values['C2']
  mean1 <- fitted.normal.values['mean1']
  mean2 <- fitted.normal.values['mean2']
  sigma1 <- fitted.normal.values['sigma1']
  sigma2 <- fitted.normal.values['sigma2']

  fun1 <- function(x) C1 * exp(-(x-mean1)**2/(2 * sigma1**2))
  fun2 <- function(x) C2 * exp(-(x-mean2)**2/(2 * sigma2**2))  
  fun3 <- function(x) C2 * exp(-(x-mean2)**2/(2 * sigma2**2)) + C1 * exp(-(x-mean1)**2/(2 * sigma1**2))

  b.granule.area <- C1 * sigma1
  a.granule.area <- C2 * sigma2
  
  a.granule.fract <- (a.granule.area / (a.granule.area + b.granule.area))
  
  parameters <- list(
    mean1=mean1,
    mean2=mean2
  )
  
  list(
    c1.fit.func=fun1,
    c2.fit.func=fun2,
    total.fit.func=fun3,
    p1=1-a.granule.fract,
    p2=a.granule.fract,
    parameters=parameters
  )  
}

find.lognormal.fit <- function(data) {
  fit <- nls(percvol.density ~ (
    C1 * dlnorm(data$bin.diameter, mean1, sigma1) + 
      C2 * dlnorm(data$bin.diameter, mean2, sigma2)
    ),
    data=data, 
    start=list(C1=30, mean1=1.5, sigma1=0.5, C2=150, mean2=3.3, sigma2=0.3),
  )
  fitted.values <- coef(fit)
  
  C1 <- fitted.values['C1']
  C2 <- fitted.values['C2']
  mean1 <- fitted.values['mean1']
  mean2 <- fitted.values['mean2']
  sigma1 <- fitted.values['sigma1']
  sigma2 <- fitted.values['sigma2']
  
  fun1 <- function(x) C1 * dlnorm(x, meanlog=mean1, sigma1)
  fun2 <- function(x) C2 * dlnorm(x, meanlog=mean2, sigma2)
  fun3 <- function(x) fun1(x) + fun2(x) 
  
  parameters <- list(
    mean1=exp(mean1),
    mean2=exp(mean2),
    sigma1=sigma1,
    sigma2=sigma2,
    p1=C1/(C1+C2),
    p2=C2/(C1+C2)
  )
  
  list(
    parameters=parameters,
    c1.fit.func=fun1,
    c2.fit.func=fun2,
    total.fit.func=fun3
  )
}


get.fitting.parameters <- function(data) {
  mixpar <- data.frame(pi=c(0.2, 0.8), mu=c(5, 20), sigma=c(2, 5))
  mixcon <- mixconstr(conpi="NONE", conmu="NONE", consigma="NONE")
  norm.fitting.results <- mix(data[,c("bin.diameter", "perc.vol")], mixpar, constr=mixcon)
  log.fitting.results <- mix(data[,c("bin.diameter", "perc.vol")], mixpar, constr=mixcon, dist="lnorm")
    
  pi1 = norm.fitting.results$parameters$pi[1]
  pi2 = norm.fitting.results$parameters$pi[2]
  
  list(
    norm.pi1=norm.fitting.results$parameters$pi[1],
    norm.pi2=norm.fitting.results$parameters$pi[2],
    norm.chisq=norm.fitting.results$chisq,
    norm.mu1=norm.fitting.results$parameters$mu[1],
    norm.mu2=norm.fitting.results$parameters$mu[2],
    norm.sigma1=norm.fitting.results$parameters$sigma[1],
    norm.sigma2=norm.fitting.results$parameters$sigma[2],
    lnorm.pi1=log.fitting.results$parameters$pi[1],
    lnorm.pi2=log.fitting.results$parameters$pi[2],
    lnorm.chisq=log.fitting.results$chisq,
    lnorm.mu1=log.fitting.results$parameters$mu[1],
    lnorm.mu2=log.fitting.results$parameters$mu[2],
    lnorm.sigma1=log.fitting.results$parameters$sigma[1],
    lnorm.sigma2=log.fitting.results$parameters$sigma[2]
  )
}

create.comparison.plot <- function(data) {
  
  ln.fit <- find.lognormal.fit(data)
  norm.fit <- find.normal.fit(data)
  
  combined.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=ln.fit$total.fit.func, color='yellow') +
    stat_function(fun=norm.fit$total.fit.func, color='purple')
  
  ypos <- max(data$percvol.density) - 0.1
  label <- sprintf("red=%.2f%%, green=%.2f%%", 100*ln.fit$p1, 100*ln.fit$p2)
  ln.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=ln.fit$total.fit.func, color='blue') +
    stat_function(fun=ln.fit$c1.fit.func, color='red') +
    stat_function(fun=ln.fit$c2.fit.func, color='green') +
    annotate("text", label=label, x=30, y=ypos)
  
  label <- sprintf("red=%.2f%%, green=%.2f%%", 100*norm.fit$p1, 100*norm.fit$p2)
  norm.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=norm.fit$total.fit.func, color='blue') +
    stat_function(fun=norm.fit$c1.fit.func, color='red') +
    stat_function(fun=norm.fit$c2.fit.func, color='green') +
    annotate("text", label=label, x=30, y=ypos)
  
  g <- arrangeGrob(ln.plot, norm.plot, combined.plot, top="test", nrow=3)
  
  g
}

fit.and.compare <- function(data, title) {

  params <- get.fitting.parameters(data)

  mu1 = params$lnorm.mu1
  s1 = params$lnorm.sigma1
  mu2 = params$lnorm.mu2
  s2 = params$lnorm.sigma2

  sd1 = sqrt(log(1 + (s1 ** 2) / (exp(2) * mu1)))
  sd2 = sqrt(log(1 + (s2 ** 2) / (exp(2) * mu2)))
  m1 = log(mu1) - (sd1 ** 2) / 2
  m2 = log(mu2) - (sd2 ** 2) / 2

  fit.ln <- nls(percvol.density ~ (
    C1 * dlnorm(data$bin.diameter, mean1, sigma1) + 
      C2 * dlnorm(data$bin.diameter, mean2, sigma2)
    ),
    data=data, 
    start=list(C1=100 * params$lnorm.pi1, mean1=m1, sigma1=sd1, C2=100*params$lnorm.pi2, mean2=m2, sigma2=sd2),
    control = list(maxiter = 5000)
  )
  fitted.values <- coef(fit.ln)

  C1.lnorm <- fitted.values['C1']
  C2.lnorm <- fitted.values['C2']
  mean1.lnorm <- fitted.values['mean1']
  mean2.lnorm <- fitted.values['mean2']
  sigma1.lnorm <- fitted.values['sigma1']
  sigma2.lnorm <- fitted.values['sigma2']
  fun1.ln <- function(x) C1.lnorm * dlnorm(x, meanlog=mean1.lnorm, sigma1.lnorm)
  fun2.ln <- function(x) C2.lnorm * dlnorm(x, meanlog=mean2.lnorm, sigma2.lnorm)
  fun3.ln <- function(x) fun1.ln(x) + fun2.ln(x) 
  lnorm.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=fun1.ln, color='blue') +
    stat_function(fun=fun2.ln, color='red') +
    stat_function(fun=fun3.ln, color='green') +
    ggtitle("Lognorm fit")


  m1 = params$norm.mu1
  m2 = params$norm.mu2
  s1 = params$norm.sigma1
  s2 = params$norm.sigma2
  fit.norm <- nls(percvol.density ~ (
    C1 * dnorm(data$bin.diameter, mean1, sigma1) + 
      C2 * dnorm(data$bin.diameter, mean2, sigma2)
    ),
    data=data, 
    start=list(C1=100 * params$norm.pi1, mean1=m1, sigma1=s1, C2=100*params$norm.pi2, mean2=m2, sigma2=s2),
  )
  fitted.values <- coef(fit.norm)
  C1 <- fitted.values['C1']
  C2 <- fitted.values['C2']
  mean1 <- fitted.values['mean1']
  mean2 <- fitted.values['mean2']
  sigma1 <- fitted.values['sigma1']
  sigma2 <- fitted.values['sigma2']
  fun1.n <- function(x) C1 * dnorm(x, mean1, sigma1)
  fun2.n <- function(x) C2 * dnorm(x, mean2, sigma2)
  fun3.n <- function(x) fun1.n(x) + fun2.n(x) 
  norm.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=fun1.n, color='blue') +
    stat_function(fun=fun2.n, color='red') +
    stat_function(fun=fun3.n, color='green') +
    ggtitle("Normal fit")


  combined.plot <- ggplot(data, aes(x=bin.diameter, y=percvol.density)) +
    geom_line() +
    stat_function(fun=fun3.ln, color='yellow') +
    stat_function(fun=fun3.n, color='purple') + 
    ggtitle("Comparison plot")

  grob <- arrangeGrob(lnorm.plot, norm.plot, combined.plot, top=title, nrow=3)

  g <- grid.arrange(grob, ncol=1)

  g
}
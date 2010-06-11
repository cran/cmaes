test.regr_stopfitness1 <- function() {
  f <- function(x, ...)
    crossprod(x)

  n <- 2
  start <- c(0, 0)
  res <- cma_es(start, f,
                lower=rep(-10, n), upper=rep(10, n),
                control=list(stopfitness=1e-5, maxit=400))

  checkTrue(res$convergence == 0)
  checkTrue(res$value < 1e-5)
}

test.regr_stopfitness2 <- function() {
  f <- function(x, ...)
    -crossprod(x)

  n <- 2
  start <- c(5, 5)
  res <- cma_es(start, f,
                lower=rep(-10, n), upper=rep(10, n),
                control=list(stopfitness=-1e-5, fnscale=-1, maxit=400))

  checkTrue(res$convergence == 0)
  checkTrue(res$value > -1e-5)
}

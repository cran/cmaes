##
## cmaES - covariance matrix adapting evolutionary strategy
##

##'
##' Global optimization procedure using a covariance matrix adapting
##' evolutionary strategy.
##'
##' Note that arguments after \code{\dots} must be matched exactly.
##' By default this function performs minimization, but it will
##' maximize if \code{control$fnscale} is negative. It can usually be
##' used as a drop in replacement for \code{optim}, but do note, that
##' no sophisticated convergence detection is included. Therefore you
##' need to Choose \code{maxit} appropriately.
##'
##' The \code{control} argument is a list that can supply any of the
##' following components:
##' \describe{
##'   \item{\code{fnscale}}{An overall scaling to be applied to the value
##'     of \code{fn} during optimization. If negative,
##'     turns the problem into a maximization problem. Optimization is
##'     performed on \code{fn(par)/fnscale}.}
##'   \item{\code{maxit}}{The maximum number of iterations. Defaults to
##'     \eqn{100*D^2}, where \eqn{D} is the dimension of the parameter space.}
##'   \item{\code{stopfitness}}{Stop if function value is smaller than or
##'     equal to \code{stopfitness}. This is the only way for the CMA-ES
##'     to \dQuote{converge}.}
##'   \item{keep.best}{return the best overall solution and not the best
##'     solution in the last population. Defaults to true.}
##'   \item{\code{sigma}}{Inital variance estimates. Can be a single
##'     number or a vector of length \eqn{D}, where \eqn{D} is the dimension
##'     of the parameter space.}
##'   \item{\code{mu}}{Population size.}
##'   \item{\code{lambda}}{Number of offspring. Must be greater than or
##'     equal to \code{mu}.}
##'   \item{\code{weights}}{Recombination weights}
##'   \item{\code{damps}}{Damping for step-size}
##'   \item{\code{cs}}{Cumulation constant for step-size}
##'   \item{\code{ccum}}{Cumulation constant for covariance matrix}
##'   \item{\code{ccov.1}}{Learning rate for rank-one update}
##'   \item{\code{ccov.mu}}{Learning rate for rank-mu update}
##'   \item{\code{diag.sigma}}{Save current step size \eqn{\sigma}{sigma}
##'     in each iteration.}
##'   \item{\code{diag.eigen}}{Save current principle components
##'     of the covariance matrix \eqn{C}{C} in each iteration.}
##'   \item{\code{diag.pop}}{Save current population in each iteration.}}
##'
##' @param par Initial values for the parameters to be optimized over.
##' @param fn A function to be minimized (or maximized), with first
##'   argument the vector of parameters over which minimization is to take
##'   place. It should return a scalar result.
##' @param \dots further arguments to be passed to \code{fn}.
##' @param lower,upper Bounds on the variables.
##' @param control A list of control parameters. See \sQuote{Details}.
##'
##' @return A list with components:
##'   \item{par}{The best set of parameters found.}
##'   \item{value}{The value of \code{fn} corresponding to \code{par}.}
##'   \item{counts}{A two-element integer vector giving the number of calls
##'     to \code{fn}. The second element is always zero for call
##'     compatibility with \code{optim}.}
##'   \item{convergence}{An integer code. \code{0} indicates successful
##'     convergence. Possible error codes are \describe{
##'       \item{\code{1}}{indicates that the iteration limit \code{maxit}
##'         had been reached.}}}
##'   \item{message}{Always set to \code{NULL}, provided for call
##'     compatibility with \code{optim}.}
##'   \item{diag}{List containing diagnostic information. Possible elements
##'     are: \describe{
##'       \item{sigma}{Vector containing the step size \eqn{\sigma}{sigma}
##'         for each iteration.}
##'       \item{eigen}{\eqn{d \times niter}{d * niter} matrix containing the
##'         principle components of the covariance matrix \eqn{C}{C}.}
##'       \item{pop}{An
##'         \eqn{d\times\mu\times niter}{d * mu * niter} array containing all
##'         populations. The last dimension is the iteration and the second
##'         dimension the individual.}}
##'    These are only present if the respective diagnostic control variable is
##'    set to \code{TRUE}.}
##'
##' @source The code is based on \file{purecmaes.m} by N. Hansen.
##'
##' @references
##' Hansen, N. (2006). The CMA Evolution Strategy: A Comparing Review. In
##'   J.A. Lozano, P. Larranga, I. Inza and E. Bengoetxea (eds.). Towards a
##'   new evolutionary computation. Advances in estimation of distribution
##'   algorithms. pp. 75-102, Springer;
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de} and
##'   David Arnu \email{david.arnu@@tu-dortmun.de}
##'
##' @title Covariance matrix adapting evolutionary strategy
##' @keywords optimize
##' @export
cma_es <- function(par, fn, ..., lower, upper, control=list()) {
  norm <- function(x)
    drop(sqrt(crossprod(x)))
  
  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }

  ## Inital solution:
  xmean <- par
  N <- length(xmean)
  ## Box constraints:
  if (missing(lower))
    lower <- rep(-Inf, N)
  else if (length(lower) == 1)  
    lower <- rep(lower, N)

  if (missing(upper))
    upper <- rep(Inf, N)
  else if (length(upper) == 1)  
    upper <- rep(upper, N)

  ## Parameters:
  trace       <- controlParam("trace", FALSE)
  fnscale     <- controlParam("fnscale", 1)
  stopfitness <- controlParam("stopfitness", -Inf)
  maxiter     <- controlParam("maxit", 100 * N^2)
  sigma       <- controlParam("sigma", 0.5)
  keep.best   <- controlParam("keep.best", TRUE)
  ## Logging options:
  diag.sigma  <- controlParam("diag.sigma", FALSE)
  diag.eigen  <- controlParam("diag.eigen", FALSE)
  diag.value   <- controlParam("diag.value", FALSE)
  diag.pop    <- controlParam("diag.pop", FALSE)
  
  ## Strategy parameter setting (defaults as recommended by Nicolas Hansen):
  lambda      <- controlParam("lambda", 4+floor(3*log(N)))
  mu          <- controlParam("mu", floor(lambda/2))
  weights     <- controlParam("weights", log(mu+1) - log(1:mu))
  weights     <- weights/sum(weights)
  mueff       <- controlParam("mueff", sum(weights)^2/sum(weights^2))
  cc          <- controlParam("ccum", 4/(N+4))
  cs          <- controlParam("cs", (mueff+2)/(N+mueff+3))
  mucov       <- controlParam("ccov.mu", mueff)
  ccov        <- controlParam("ccov.1",
                              (1/mucov) * 2/(N+1.4)^2
                              + (1-1/mucov) * ((2*mucov-1)/((N+2)^2+2*mucov)))
  damps       <- controlParam("damps",
                              1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs)

  ## Safety checks:
  stopifnot(length(upper) == N)  
  stopifnot(length(lower) == N)
  stopifnot(all(lower < upper))
  stopifnot(length(sigma) == 1)

  ## Bookkeeping variables for the best solution found so far:
  best.fit <- Inf
  best.par <- NULL

  ## Preallocate logging structures:
  if (diag.sigma)
    sigma.log <- numeric(maxiter)
  if (diag.eigen)
    eigen.log <- matrix(0, nrow=maxiter, ncol=N)
  if (diag.value)
    value.log <- numeric(maxiter)
  if (diag.pop)
    pop.log <- array(0, c(N, mu, maxiter))
  
  ## Initialize dynamic (internal) strategy parameters and constants
  pc <- rep(0.0, N)
  ps <- rep(0.0, N)
  B <- diag(N)
  D <- diag(N)
  BD <- B %*% D
  C <- BD %*% t(BD)

  chiN <- sqrt(N) * (1-1/(4*N)+1/(21*N^2))
  
  iter <- 0L      ## Number of iterations
  counteval <- 0L ## Number of function evaluations
  cviol <- 0L     ## Number of constraint violations
  msg <- NULL     ## Reason for terminating
  
  ## Preallocate work arrays:
  arx <- matrix(0.0, nrow=N, ncol=lambda)
  arz <- matrix(0, nrow=N, ncol=lambda)
  arfitness <- numeric(lambda)
  while (iter < maxiter) {
    iter <- iter + 1L

    if (diag.sigma)
      sigma.log[iter] <- sigma

    ## Generate new population:
    arz <- matrix(rnorm(N*lambda), ncol=lambda)
    for (k in 1:lambda) {
      x <- xmean + sigma * (BD %*% arz[,k])
      ## Calculate 'valid' x inside box constraints and penalty:
      vx <- pmin(pmax(x, lower), upper)
      pen <- 1 + norm(x - vx)        
      if (pen > 1.0) {
        cviol <- cviol + 1
        if (pen == Inf)
          pen <- 1e200
      }
      ## Calculate fitness:
      arfitness[k] <- fn(as.vector(vx), ...) * pen * fnscale
      counteval <- counteval + 1L;
      ## Save point:
      arx[,k] <- x
    }
    
    ## Order fitness:
    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]

    aripop <- arindex[1:mu]
    selx <- arx[,aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[,aripop]
    zmean <- drop(selz %*% weights)

    ## Save selected x value:
    if (diag.pop) pop.log[,,iter] <- selx

    ## Cumulation: Update evolutionary paths
    ps <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN) < (1.4 + 2/(N+1)))
    pc <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * drop(BD %*% zmean)

    ## Adapt Covariance Matrix:
    BDz <- BD %*% selz
    C <- (1-ccov) * C + ccov * (1/mucov) *
      (pc %o% pc + (1-hsig) * cc*(2-cc) * C) +
        ccov * (1-1/mucov) * BDz %*% diag(weights) %*% t(BDz)
    
    ## Adapt step size sigma:
    sigma <- sigma * exp((norm(ps)/chiN - 1)*cs/damps)
    
    e <- eigen(C, symmetric=TRUE)
    if (diag.eigen)
      eigen.log[iter,] <- rev(sort(e$values))

    if (!all(e$values >= sqrt(.Machine$double.eps) * abs(e$values[1]))) {      
      msg <- "Covariance matrix 'C' is numerically not positive definite."
      break
    }

    B <- e$vectors
    D <- if (length(e$values) > 1L)
      diag(sqrt(e$values))
    else
      as.matrix(sqrt(e$values))
    BD <- B %*% D

    if (keep.best && arfitness[1] < best.fit) {
      best.fit <- arfitness[1] 
      best.par <- arx[, arindex[1]]
    }

    ## break if fit:
    if (arfitness[1] <= stopfitness * fnscale) {
      msg <- "Stop fitness reached."
      break
    }

    if (diag.value)
      value.log[iter] <- arfitness[1]
    
    ## Escape from flat-land:
    if (arfitness[1] == arfitness[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) { 
      sigma <- sigma * exp(0.2+cs/damps); 
      warning("Flat fitness function. Increasing sigma.")
    }
    if (trace)
      message(sprintf("Iteration %i of %i: current fitness %f",
                      iter, maxiter, arfitness[1] * fnscale))
  }
  cnt <- vector("integer", 2L)
  cnt[1] <- as.integer(counteval)
  cnt[2] <- NA
  names(cnt) <- c("function", "gradient")

  ## Currently only the last 'best' solution is returned, this may not
  ## be the best one found over all generations.
  if (!keep.best) {
    best.fit <- arfitness[1]
    best.par <- arx[, arindex[1]]
  }

  diag <- list()
  ## Subset diagnostic data to only include those iterations which
  ## where actually performed.
  if (diag.value) diag$value <- value.log[1:iter]
  if (diag.sigma) diag$sigma <- sigma.log[1:iter]
  if (diag.eigen) diag$eigen <- eigen.log[1:iter,]
  if (diag.pop)   diag$pop   <- pop.log[,,1:iter]
  
  res <- list(par=best.par,
              value=best.fit * fnscale,
              counts=cnt,
              convergence=ifelse(iter >= maxiter, 1L, 0L),
              message=msg,
              constr.violations=cviol,
              diag=diag
              )
  nm <- names(par)
  if (!is.null(nm)) 
    names(res$par) <- nm  
  class(res) <- "cma_es.result"
  return(res)
}

##' @rdname cma_es
##' @export
cmaES <- function(...) {
  .Deprecated("cma_es")
  cma_es(...)
}

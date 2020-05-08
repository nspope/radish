radish_optimize <- function(f, g, s, S, theta = rep(0, ncol(s$x)), leverage = TRUE, nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6))
{
  .fcall <<- 0 #debug
  problem <- BoxConstrainedNewton(theta, 
                           function(par, gradient, hessian) 
                           {
                             .fcall <<- .fcall + 1 #debug
                               radish_algorithm(f = f, g = g, s = s, S = S, theta = c(par), 
                                                gradient = gradient, 
                                                hessian = hessian, 
                                                partial = FALSE, 
                                                nonnegative = nonnegative)
                           },
                           control = control)
  theta <- problem$par
  fit   <- radish_algorithm(f = f, g = g, s = s, S = S, theta = theta,
                            gradient = TRUE, hessian = TRUE, partial = leverage,
                            nonnegative = nonnegative)

  # calculate leverage for genetic distance and spatial covariates
  if (leverage)
  {
    ihess      <- MASS::ginv(fit$hessian) 
    leverage_S <- -matrix(fit$partial_S, length(S), ncol(s$x)) %*% ihess
    leverage_S[upper.tri(S),] <- 0
    leverage_S <- array(leverage_S, dim = c(nrow(S), ncol(S), length(theta)))
    leverage_X <- array(NA, dim = dim(fit$partial_X))
    for (k in 1:ncol(s$x))
      leverage_X[,,k] <- -fit$partial_X[,,k] %*% ihess
  }

  if (validate)
  {
    silence <- function(control) { control$verbose = FALSE; control }

    num_leverage_S <- 
      array(t(numDeriv::jacobian(function(S) 
                                 radish_main(f = f,
                                             g = g, 
                                             s = s,
                                             S = S, 
                                             theta = theta,
                                             nonnegative = nonnegative,
                                             control = silence(control))$theta, 
                                 S, method = "simple")), 
            dim = dim(leverage_S))

    num_leverage_X <- array(NA, dim = dim(leverage_X))
    for (k in 1:ncol(s$x))
      num_leverage_X[,,k] <- 
        t(numDeriv::jacobian(function(z) {
                               x[,k] <- z 
                               radish_main(f = f,
                                           g = g, 
                                           s = s,
                                           S = S, 
                                           theta = theta,
                                           nonnegative = nonnegative,
                                           control = silence(control))$theta }, 
                             x[,k], method = "simple"))
  }

  if (fit$boundary)
    warning("Optimum for subproblem is on boundary (e.g. no spatial genetic structure): cannot optimize theta. Try different starting values.")
  else
  {
    ztable <- matrix(0, length(theta), 4)
    colnames(ztable) <- c("Est", "StdErr", "Z", "pval")
    rownames(ztable) <- s$covariates
    ztable[,"Est"]   <- theta
    ztable[,"StdErr"] <- diag(solve(fit$hessian))
    ztable[,"Z"]      <- ztable[,"Est"]/ztable[,"StdErr"]
    ztable[,"pval"]   <- pmin(2*(1 - pnorm(abs(ztable[,"Z"]))), 1)
  }

  list(fcall          = .fcall, #DEBUG
       fit            = fit,
       theta          = theta,
       ztable         = ztable,
       AIC            = 2*fit$objective + 2*length(theta) + 2*length(fit$phi),
       phi            = fit$phi,
       loglikelihood  = -fit$objective,
       gradient       = -fit$gradient,
       hessian        = -fit$hessian,
       leverage_S     = if(!leverage) NULL else leverage_S,
       leverage_X     = if(!leverage) NULL else leverage_X,
       num_leverage_S = if(!validate) NULL else num_leverage_S,
       num_leverage_X = if(!validate) NULL else num_leverage_X)
}

radish_optimize2 <- function(f, g, s, S, theta = rep(0, ncol(s$x)), leverage = TRUE, nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6))
{
  .fcall <<- 0 #debug
  problem <- BFGS(theta, 
                           function(par, gradient, hessian) 
                           {
                             .fcall <<- .fcall + 1 #debug
                               radish_algorithm(f = f, g = g, s = s, S = S, theta = c(par), 
                                                gradient = gradient, 
                                                hessian = hessian, 
                                                partial = FALSE, 
                                                nonnegative = nonnegative)
                           },
                           control = control)
  theta <- problem$par
  fit   <- radish_algorithm(f = f, g = g, s = s, S = S, theta = theta,
                            gradient = TRUE, hessian = TRUE, partial = leverage,
                            nonnegative = nonnegative)

  # calculate leverage for genetic distance and spatial covariates
  if (leverage)
  {
    ihess      <- MASS::ginv(fit$hessian) 
    leverage_S <- -matrix(fit$partial_S, length(S), ncol(s$x)) %*% ihess
    leverage_S[upper.tri(S),] <- 0
    leverage_S <- array(leverage_S, dim = c(nrow(S), ncol(S), length(theta)))
    leverage_X <- array(NA, dim = dim(fit$partial_X))
    for (k in 1:ncol(s$x))
      leverage_X[,,k] <- -fit$partial_X[,,k] %*% ihess
  }

  if (validate)
  {
    silence <- function(control) { control$verbose = FALSE; control }

    num_leverage_S <- 
      array(t(numDeriv::jacobian(function(S) 
                                 radish_main(f = f,
                                             g = g, 
                                             s = s,
                                             S = S, 
                                             theta = theta,
                                             nonnegative = nonnegative,
                                             control = silence(control))$theta, 
                                 S, method = "simple")), 
            dim = dim(leverage_S))

    num_leverage_X <- array(NA, dim = dim(leverage_X))
    for (k in 1:ncol(s$x))
      num_leverage_X[,,k] <- 
        t(numDeriv::jacobian(function(z) {
                               x[,k] <- z 
                               radish_main(f = f,
                                           g = g, 
                                           s = s,
                                           S = S, 
                                           theta = theta,
                                           nonnegative = nonnegative,
                                           control = silence(control))$theta }, 
                             x[,k], method = "simple"))
  }

  if (fit$boundary)
    warning("Optimum for subproblem is on boundary (e.g. no spatial genetic structure): cannot optimize theta. Try different starting values.")

  list(fcall          = .fcall, #DEBUG
       fit            = fit,
       theta          = theta,
       phi            = fit$phi,
       loglikelihood  = fit$objective,
       gradient       = fit$gradient,
       hessian        = fit$hessian,
       leverage_S     = if(!leverage) NULL else leverage_S,
       leverage_X     = if(!leverage) NULL else leverage_X,
       num_leverage_S = if(!validate) NULL else num_leverage_S,
       num_leverage_X = if(!validate) NULL else num_leverage_X)
}

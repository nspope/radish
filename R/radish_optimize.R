setRefClass("FunctionCall", fields = list(count = "integer"))

#' Optimize a parameterized conductance surface
#'
#' Calculates likelihood, gradient, hessian, and partial derivatives of a
#' parameterized conductance surface, given a function mapping spatial data to
#' conductance and a function mapping resistance distance (covariance) to
#' genetic distance.
#'
#' @param f A function of class 'conductance_model'
#' @param g A function of class 'measurement_model'
#' @param s An object of class 'radish_graph'
#' @param S A matrix of observed genetic distances
#' @param nu Number of genetic markers (potentially used by 'g')
#' @param theta Starting values for optimization
#' @param leverage Compute influence measures and leverage?
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#' @param optimizer The optimization algorithm to use: "newton" uses the exact Hessian, so computational cost grows linearly with the number of parameters; while "bfgs" uses an approximation with much reduced cost (but slower overall convergence)
#' @param control A list containing options for the optimization routine (see ?NewtonRaphsonControl for list)
#'
#' @return An object of class 'radish'
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' fit_nnls <- radish(radish::loglinear_conductance, radish::leastsquares, surface, melip.Fst)
#' summary(fit_nnls)
#'
#' @export

radish <- function(f, g, s, S, nu = NULL, theta = rep(0, ncol(s$x)), leverage = TRUE, nonnegative = TRUE, validate = FALSE, optimizer = c("newton", "bfgs"), control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6))
{
  optimizer <- match.arg(optimizer)
  fcalls    <- new("FunctionCall", count = 0L)
  optfn     <- function(par, gradient, hessian)
  {
    fcalls$count <- fcalls$count + 1L
    radish_algorithm(f = f, g = g, s = s, S = S, nu = nu, theta = c(par), 
                     gradient = gradient, 
                     hessian = hessian, 
                     partial = FALSE, 
                     nonnegative = nonnegative)
  }

  if (optimizer == "newton")
    problem <- BoxConstrainedNewton(theta, 
                                    optfn,
                                    control = control)
  else if (optimizer == "bfgs")
    problem <- BoxConstrainedBFGS(theta, 
                                  optfn,
                                  control = control)

  theta <- problem$par
  fit   <- radish_algorithm(f = f, g = g, s = s, S = S, nu = nu, theta = theta,
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
    colnames(ztable)      <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(ztable)      <- s$covariates
    ztable[,"Estimate"]   <- theta
    ztable[,"Std. Error"] <- sqrt(diag(solve(fit$hessian)))
    ztable[,"z value"]    <- ztable[,"Estimate"]/ztable[,"Std. Error"]
    ztable[,"Pr(>|z|)"]   <- pmin(2*(1 - pnorm(abs(ztable[,"z value"]))), 1)

    vcor <- cov2cor(solve(fit$hessian))
    rownames(vcor) <- colnames(vcor) <- s$covariates
  }

  out <- list(call           = match.call(),
              response       = S,
              fcalls         = fcalls$count,
              iters          = problem$iters,
              boundary       = fit$boundary,
              fit            = fit,
              theta          = theta,
              phi            = fit$phi[,1],
              ztable         = ztable,
              vcor           = as.dist(vcor),
              aic            = 2*fit$objective + 2*length(theta) + 2*length(fit$phi),
              df             = (!fit$boundary) * length(theta) + length(fit$phi),
              loglikelihood  = -fit$objective,
              gradient       = -fit$gradient,
              hessian        = -fit$hessian,
              leverage_S     = if(!leverage) NULL else leverage_S,
              leverage_X     = if(!leverage) NULL else leverage_X,
              num_leverage_S = if(!validate) NULL else num_leverage_S,
              num_leverage_X = if(!validate) NULL else num_leverage_X)
  class(out) <- "radish"
  out
}

print.radish <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Conductance surface estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (nrow(x$ztable) && !x$boundary)
  {
    cat("Coefficients:\n")
    print.default(format(x$ztable[,"Estimate"], digits = digits), print.gap = 2L, quote = FALSE)
  }
  else if (x$boundary)
  {
    cat("Model fit is on boundary (e.g. no genetic structure), no coefficients\n")
  }
  else
  {
    cat("No coefficients\n")
  }
  cat("\n")
  cat("Loglikelihood:", x$loglikelihood, "   Df:", x$df, "   AIC:", x$aic, "\n")
}

summary.radish <- function(x, ...)
{
  out <- list(boundary = x$boundary,
              ztable   = x$ztable,
              vcor     = x$vcor,
              loglikelihood = x$loglikelihood,
              aic      = x$aic,
              fcalls   = x$fcalls,
              iters    = x$iters,
              df       = x$df,
              gradnorm = sqrt(c(x$gradient %*% x$gradient)),
              call     = x$call,
              phi      = x$phi)
  class(out) <- "summary.radish"
  out
}

print.summary.radish <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...)
{
  cat("Conductance surface estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Loglikelihood:", x$loglikelihood, "   Df:", x$df,  "\n")
  cat("AIC:", x$aic, "\n\n")
  cat("Number of function calls:", x$fcalls, "\n")
  cat("Number of Newton-Raphson steps:", x$iters, "\n")
  cat("Norm of gradient at MLE:", x$gradnorm, "\n\n")
  if (length(x$phi))
  {
    cat("Nuisance parameters:\n")
    print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  if (nrow(x$ztable) && !x$boundary)
  {
    cat("Coefficients:\n")
    printCoefmat(x$ztable, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    if (nrow(x$ztable) > 1)
    {
      cat("\n")
      cat("Correlation of Coefficients:\n")
      print(x$vcor)
    }
  }
  else if (x$boundary)
  {
    cat("Model fit is on boundary (e.g. no genetic structure), no valid coefficients\n")
  }
  else
  {
    cat("No coefficients\n")
  }
}

fitted.radish <- function(x, type = c("response", "distance", "covariance"), ...)
{
  type <- match.arg(type)
  if (type == "response")
    as.matrix(x$fit$fitted)
  else if (type == "distance")
    dist_from_cov(as.matrix(x$fit$covariance))
  else if (type == "covariance")
    as.matrix(x$fit$covariance)
}

simulate.radish <- function(x, nsim = 1, method = c("permutation", "parametric"), ...)
{
  method <- match.arg(method)
  if (method == "parametric")
  {
    stop("Parametric simulation not yet supported")
  } 
  else if (method == "permutation") 
  {
    fit    <- fitted(x)
    resid  <- x$response - fit
    sims   <- array(NA, c(nrow(resid), ncol(resid), nsim))
    for (i in 1:nsim)
    {
      ind <- sample(1:nrow(resid))
      sims[,,i] <- fit + resid[ind,ind]
    }
  }
  if (nsim == 1) sims[,,1] else sims
}

# at some point get a non-finicky version of BFGS
# needs debugging
radish2 <- function(f, g, s, S, nu = NULL, theta = rep(0, ncol(s$x)), leverage = TRUE, nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6))
{
  fcalls  <- new("FunctionCall", count = 0L)
  problem <- BoxConstrainedBFGS(theta, 
                                function(par, gradient, hessian) 
                                {
                                  fcalls$count <- fcalls$count + 1L
                                  radish_algorithm(f = f, g = g, s = s, S = S, nu = nu, theta = c(par), 
                                                   gradient = gradient, 
                                                   hessian = hessian, 
                                                   partial = FALSE, 
                                                   nonnegative = nonnegative)
                                },
                                control = control)
  theta <- problem$par
  fit   <- radish_algorithm(f = f, g = g, s = s, S = S, nu = nu, theta = theta,
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
    colnames(ztable)      <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(ztable)      <- s$covariates
    ztable[,"Estimate"]   <- theta
    ztable[,"Std. Error"] <- diag(solve(fit$hessian))
    ztable[,"z value"]    <- ztable[,"Estimate"]/ztable[,"Std. Error"]
    ztable[,"Pr(>|z|)"]   <- pmin(2*(1 - pnorm(abs(ztable[,"z value"]))), 1)

    vcor <- cov2cor(solve(fit$hessian))
    rownames(vcor) <- colnames(vcor) <- s$covariates
  }

  out <- list(call           = match.call(),
              response       = S,
              fcalls         = fcalls$count,
              iters          = problem$iters,
              boundary       = fit$boundary,
              fit            = fit,
              theta          = theta,
              phi            = fit$phi[,1],
              ztable         = ztable,
              vcor           = as.dist(vcor),
              aic            = 2*fit$objective + 2*length(theta) + 2*length(fit$phi),
              df             = (!fit$boundary) * length(theta) + length(fit$phi),
              loglikelihood  = -fit$objective,
              gradient       = -fit$gradient,
              hessian        = -fit$hessian,
              leverage_S     = if(!leverage) NULL else leverage_S,
              leverage_X     = if(!leverage) NULL else leverage_X,
              num_leverage_S = if(!validate) NULL else num_leverage_S,
              num_leverage_X = if(!validate) NULL else num_leverage_X)
  class(out) <- "radish"
  out
}

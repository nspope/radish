setRefClass("FunctionCall", fields = list(count = "integer"))

#' Optimize a parameterized conductance surface
#'
#' Uses maximum likelihood to fit a parameterized conductance surface to genetic data,
#' given a function mapping spatial data to conductance (a "conductance model")
#' and a function mapping resistance distance (covariance) to genetic distance
#' (a "measurement model").
#'
#' @param formula A formula with the name of a matrix of observed genetic distances on the lhs, and covariates in the creation of \code{data} on the rhs
#' @param data An object of class \code{radish_graph} (see \code{\link{conductance_surface}})
#' @param conductance_model A function of class \code{radish_conductance_model_factory} (see \code{\link{radish_conductance_model_factory}})
#' @param measurement_model A function of class \code{radish_measurement_model} (see \code{\link{radish_measurement_model}})
#' @param nu Number of genetic markers (potentially used by \code{measurement_model})
#' @param theta Starting values for optimization
#' @param leverage Compute influence measures and leverage?
#' @param nonnegative Force regression-like \code{measurement_model} to have nonnegative slope?
#' @param conductance If \code{TRUE}, edge conductance is the sum of cell conductances; otherwise edge conductance is the inverse of the sum of cell resistances (NOT USED; TODO)
#' @param optimizer The optimization algorithm to use: \code{newton} uses the exact Hessian, with computational cost that grows linearly with the number of parameters; while \code{bfgs} uses an approximation with much reduced cost (but slower overall convergence)
#' @param control A list containing options for the optimization routine (see \code{\link{NewtonRaphsonControl}} for list)
#' @param validate Numerical validation of leverage via package \code{numDeriv} (very slow, use for debugging small examples)
#'
#' @details By "parameterized conductance surface", what is meant is a model where per-vertex conductance (and thus resistance distance) is a function of spatial covariates. The choice of function is referred to in this package as the "conductance model". The inverse problem (and the purpose of this package) is to estimate the parameters of the conductance model, by relating the (unknown, modelled) resistance distance to observed genetic dissimilarity via a probability model (referred to as the "measurement model" throughout this package).
#'
#' For example, a log-linear choice of conductance model is:
#'
#'   \code{conductance[i] = exp(covariates[i,] \%*\% theta)}
#'
#' where \code{theta} are unknown parameters and \code{covariates} is a design matrix. \code{radish} estimates \code{theta} as well as any nuisance parameters associated with the measurement model.
#'
#' If the fit is on the boundary (e.g. no spatial genetic structure) or is the null model of isolation-by-distance, the resulting object will not contain influence/leverage/gradient/hessian.
#'
#' @return An object of class \code{radish}
#'
#' @examples
#'
#' library(raster)
#' 
#' data(melip)
#' 
#' # scaling spatial covariates helps avoid numeric overflow
#' covariates <- raster::stack(list(altitude = raster::scale(melip.altitude), 
#'                                  forestcover = raster::scale(melip.forestcover)))
#' 
#' # 
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' fit_nnls <- radish(melip.Fst ~ altitude * forestcover, data = surface, 
#'                    radish::loglinear_conductance, radish::leastsquares)
#' summary(fit_nnls)
#' 
#' # different "measurement_model" that incorporates dependence
#' # among pairwise measurements
#' fit_mlpe <- radish(melip.Fst ~ altitude * forestcover, data = surface, 
#'                    radish::loglinear_conductance, radish::mlpe)
#' summary(fit_mlpe)
#'
#' # conductance surface with 95% CI
#' fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)
#' 
#' # test for an interaction using a likelihood ratio test
#' fit_mlpe_interaction <- radish(melip.Fst ~ forestcover * altitude, data = surface, 
#'                                radish::loglinear_conductance, radish::mlpe)
#' anova(fit_mlpe, fit_mlpe_interaction)
#' 
#' # test against null model of IBD using a LRT
#' fit_mlpe_ibd <- radish(melip.Fst ~ 1, data = surface, 
#'                        radish::loglinear_conductance, radish::mlpe)
#' anova(fit_mlpe, fit_mlpe_ibd)
#'
#' @export

radish <- function(formula, 
                   data,
                   conductance_model = radish::loglinear_conductance, 
                   measurement_model = radish::mlpe, 
                   nu = NULL, 
                   theta = NULL,
                   leverage = TRUE, 
                   nonnegative = TRUE, 
                   conductance = TRUE, 
                   optimizer = c("newton", "bfgs"), 
                   control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6), 
                   validate = FALSE)
{
  stopifnot(class(formula) == "formula")
  stopifnot(class(data) == "radish_graph")
  stopifnot(class(conductance_model) == "radish_conductance_model_factory")
  stopifnot(class(measurement_model) == "radish_measurement_model")

  # get response, remove lhs from formula
  terms    <- terms(formula)
  vars     <- as.character(attr(terms, "variables"))[-1]
  response <- attr(terms, "response")
  S        <- if(response) get(vars[attr(terms, "response")], parent.frame())
              else stop("'formula' must have genetic distance matrix on lhs")
  is_ibd   <- length(vars) == 1
  formula  <- if (!is_ibd) reformulate(attr(terms, "term.labels"))
              else formula(~1)

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model <- conductance_model(formula, data$x) #TODO: is_ibd here and have factory modify accordingly

  # initialize theta
  default <- attr(conductance_model, "default")
  if (is.null(theta))
    theta <- default
  else
    stopifnot(length(theta) == length(default))
  names(theta) <- names(default)

  optimizer <- match.arg(optimizer)
  fcalls    <- new("FunctionCall", count = 0L)
  optfn     <- function(par, gradient, hessian)
  {
    fcalls$count <- fcalls$count + 1L
    radish_algorithm(f = conductance_model, 
                     g = measurement_model, 
                     s = data,
                     S = S, 
                     nu = nu, 
                     theta = c(par), 
                     gradient = gradient, 
                     hessian = hessian, 
                     partial = FALSE, 
                     nonnegative = nonnegative)
  }

  if (!is_ibd)
  {
    if (optimizer == "newton")
      problem <- BoxConstrainedNewton(theta, 
                                      optfn,
                                      control = control)
    else if (optimizer == "bfgs")
      problem <- BoxConstrainedBFGS(theta, 
                                    optfn,
                                    control = control)
    iters <- problem$iters
    theta <- c(problem$par)
    names(theta) <- names(default)
  }
  else
  {
    iters <- 0L
    theta <- c(0) #TODO: this will not necessarily equate to IBD, see above
  }

  fit <- radish_algorithm(f = conductance_model, g = measurement_model, 
                          s = data, S = S, nu = nu, theta = theta,
                          gradient = TRUE, hessian = TRUE, partial = leverage,
                          nonnegative = nonnegative)

  fit$response <- S

  if (fit$boundary)
    warning("Optimum for subproblem is on boundary (e.g. no spatial genetic structure): cannot optimize theta.\nTry different starting values.")
  no_coef <- fit$boundary || is_ibd 

  # calculate leverage for genetic distance and spatial covariates
  leverage <- leverage && !no_coef
  if (leverage)
  {
    ihess      <- MASS::ginv(fit$hessian) 
    leverage_S <- -matrix(fit$partial_S, length(S), length(theta)) %*% ihess
    leverage_S[upper.tri(S),] <- 0
    leverage_S <- array(leverage_S, dim = c(nrow(S), ncol(S), length(theta)))
    leverage_X <- array(NA, dim = dim(fit$partial_X))
    for (k in 1:length(theta))
      leverage_X[,,k] <- -fit$partial_X[,,k] %*% ihess
  }

  # numerical validation of derivatives
  validate <- validate && !no_coef
  if (validate)
  {
    #TODO update with new API
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
    for (k in 1:length(theta))
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

  out <- list(call           = match.call(),
              formula        = formula,
              dim            = c("vertices" = nrow(data$x),
                                 "focal"    = nrow(S),
                                 "edge"     = ncol(data$adj)),
              cost           = c("newton_steps"   = iters,
                                 "function_calls" = fcalls$count + 1),
              submodels      = list("f" = conductance_model,
                                    "g" = measurement_model),
              fit            = fit,
              loglik         = -fit$objective,
              df             = (!no_coef)*length(theta) + length(fit$phi),
              aic            = 2*fit$objective + 2*(!no_coef)*length(theta) + 2*length(fit$phi),
              mle            = list("theta"    = if(no_coef) NULL else theta,
                                    "gradient" = if(no_coef) NULL else -fit$gradient,
                                    "hessian"  = if(no_coef) NULL else -fit$hessian),
              leverage       = list("S" = if(!leverage) NULL else leverage_S,
                                    "X" = if(!leverage) NULL else leverage_X,
                                    "validate" = if(!validate || !leverage) NULL 
                                                 else list("S" = num_leverage_S,
                                                           "X" = num_leverage_X))
              )
  class(out) <- "radish"
  out
}

print.radish <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Conductance surface estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!x$fit$boundary && !is.null(x$mle$theta))
  {
    cat("Coefficients:\n")
    print.default(format(x$mle$theta, digits = digits), print.gap = 2L, quote = FALSE)
  }
  else
  {
    cat("No coefficients\n")
  }
  cat("\n")
  cat("Loglikelihood:", x$loglik, paste0("(", x$df), "degrees freedom)   AIC:", x$aic, "\n")
}

summary.radish <- function(x, ...)
{
  tol <- sqrt(.Machine$double.eps) #for checking singularity

  no_coef <- x$fit$boundary || is.null(x$mle$theta)
  if (!no_coef)
  {
    ztable <- matrix(0, length(x$mle$theta), 4)
    colnames(ztable)      <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(ztable)      <- names(x$mle$theta)
    ztable[,"Estimate"]   <- x$mle$theta
    ztable[,"Std. Error"] <- sqrt(diag(solve(x$fit$hessian)))
    ztable[,"z value"]    <- ztable[,"Estimate"]/ztable[,"Std. Error"]
    ztable[,"Pr(>|z|)"]   <- pmin(2*(1 - pnorm(abs(ztable[,"z value"]))), 1)

    ehess <- eigen(x$fit$hessian)

    if (any(ehess$values < 0))
      warning("Hessian matrix has negative eigenvalues: possibly a saddle point")
    if (any(abs(ehess$values) < tol * max(abs(ehess$values))))
      warning("Hessian matrix is singular or nearly singular: model is probably non-identifiable")

    vcov <- ehess$vectors %*% diag(1/ehess$values, nrow = nrow(x$fit$hessian)) %*% t(ehess$vectors)
    vcor <- cov2cor(vcov)
    rownames(vcor) <- colnames(vcor) <- 
      rownames(vcov) <- colnames(vcov) <- 
        rownames(ztable)
  }

  out <- list(boundary      = x$fit$boundary,
              phi           = x$fit$phi[,1],
              ztable        = if (no_coef) NULL else ztable,
              vcor          = if (no_coef) NULL else vcor,
              vcov          = if (no_coef) NULL else vcov,
              gradnorm      = if (no_coef) NA else sqrt(sum(x$mle$gradient^2)),
              loglik        = x$loglik,
              df            = x$df,
              aic           = x$aic,
              fcalls        = x$cost["function_calls"],
              iters         = x$cost["newton_steps"],
              call          = x$call,
              dim           = x$dim
              )
  class(out) <- "summary.radish"
  out
}

print.summary.radish <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...)
{
  cat("Conductance surface with", x$dim["vertices"], "vertices", 
      paste0("(", x$dim["focal"]), "focal) estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Loglikelihood:", x$loglik, paste0("(", x$df), "degrees freedom)\nAIC:", x$aic, "\n\n")
  cat("Number of function calls:", x$fcalls, "\n")
  cat("Number of Newton-Raphson steps:", x$iters, "\n")
  cat("Norm of gradient at MLE:", x$gradnorm, "\n\n")
  if (length(x$phi))
  {
    cat("Nuisance parameters:\n")
    print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  if (!x$boundary && !is.null(x$ztable))
  {
    cat("Coefficients:\n")
    printCoefmat(x$ztable, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    if (nrow(x$ztable) > 1)
    {
      cat("\n")
      cat("Correlation of Coefficients:\n")
      print(as.dist(x$vcor))
    }
  }
  else if (x$boundary)
  {
    cat("Model fit is on boundary (e.g. no genetic structure), no meaningful coefficients\n")
  }
  else
  {
    cat("No coefficients\n")
  }
}

coef.radish <- function(x)
{
  if (!is.null(x$mle$theta))
    return(x$mle$theta)
  else
    return(c())
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
    resid  <- x$fit$response - fit
    sims   <- array(NA, c(nrow(resid), ncol(resid), nsim))
    for (i in 1:nsim)
    {
      ind <- sample(1:nrow(resid))
      sims[,,i] <- fit + resid[ind,ind]
    }
  }
  if (nsim == 1) sims[,,1] else sims
}

anova.radish <- function(object, alternative)
{
  stopifnot(class(object) == "radish" && class(alternative) == "radish")
  stopifnot(!object$fit$boundary && !alternative$fit$boundary)

  if (object$df >= alternative$df)
  {
    full <- object
    reduced <- alternative
  } 
  else
  {
    full <- alternative
    reduced <- object
  }

  form_reduced <- paste("Null:", paste(reduced$formula, collapse = " "))
  form_full    <- paste("Alt:", paste(full$formula, collapse = " "))

  Chisq <- 2 * (full$loglik - reduced$loglik)
  Df    <- full$df - reduced$df
  P     <- pchisq(Chisq, Df, lower = FALSE)
  Ll    <- c(reduced$loglik, full$loglik)
  Np    <- c(reduced$df, full$df)

  out <- cbind("logLik" = Ll, 
               "Df" = Np, 
               "ChiSq" = c(NA, Chisq), 
               "Df(ChiSq)" = c(NA, Df), 
               "Pr(>Chi)" = c(NA, P))
  rownames(out) <- c("Null", "Alt")

  attr(out, "heading") <- c("Likelihood ratio test",
                           form_reduced, form_full)
  class(out) <- "anova"
  out
}

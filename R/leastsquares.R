#' (Nonnegative) least squares
#'
#' A function of class \code{measurement_model} that calculates likelihood,
#' gradient, hessian, and partial derivatives of nuisance parameters and the
#' Laplacian generalized inverse, using nonnegative least squares.
#'
#' @param E A submatrix of the generalized inverse of the graph Laplacian (e.g. a covariance matrix)
#' @param S A matrix of observed genetic distances
#' @param phi Nuisance parameters (see details)
#' @param nu Unused
#' @param gradient Compute gradient of negative loglikelihood with regard to \code{phi}?
#' @param hessian Compute Hessian matrix of negative loglikelihood with regard to \code{phi}?
#' @param partial Compute second partial derivatives of negative loglikelihood with regard to \code{phi}, \code{E}, \code{S}?
#' @param nonnegative Force slope to be nonnegative?
#' @param validate Numerical validation via package \code{numDeriv} (very slow, use for debugging small examples)
#'
#' @details The nuisance parameters \code{phi} are the intercept ("alpha"), slope ("beta"), and negative log residual
#' standard deviation ("tau") of the least squares regression. If not supplied, \code{phi} 
#' is estimated via maximum likelihood by \code{nlme::gls}.
#'
#' TODO: formula
#'
#' @seealso \code{\link{radish_measurement_model}}
#'
#' @return A list containing:
#'  \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
#'  \item{objective}{(if \code{objective}) the negative loglikelihood}
#'  \item{fitted}{((if \code{objective}) a matrix of expected genetic distances among target vertices}
#'  \item{boundary}{(if \code{objective}) is the MLE on the boundary (e.g. no genetic structure)?}
#'  \item{gradient}{(if \code{gradient}) gradient of negative loglikelihood with respect to phi}
#'  \item{hessian}{(if \code{hessian}) Hessian matrix of the negative loglikelihood with respect to phi}
#'  \item{gradient_E}{(if \code{partial}) gradient with respect to the generalized inverse of the graph Laplacian}
#'  \item{partial_E}{(if \code{partial}) Jacobian of \code{gradient_E} with respect to phi}
#'  \item{partial_S}{(if \code{partial}) Jacobian of \code{gradient} with respect to S}
#'  \item{jacobian_E}{(if \code{partial}) a function used for reverse algorithmic differentiation}
#'  \item{jacobian_S}{(if \code{partial}) a function used for reverse algorithmic differentiation}
#'
#' @examples
#'
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' # inverse of graph Laplacian at null model (IBD) 
#' laplacian_inv <- radish_distance(theta = matrix(0, 1, 2), 
#'                                  formula = ~forestcover + altitude,
#'                                  data = surface,
#'                                  radish::loglinear_conductance, 
#'                                  covariance = TRUE)$covariance[,,1]
#' 
#' leastsquares(laplacian_inv, melip.Fst) #without 'phi': return MLE of phi
#' leastsquares(laplacian_inv, melip.Fst, phi = c(0., 0.5, -0.1))
#'
#' @export

leastsquares <- function(E, S, phi, nu = NULL, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  symm <- function(X) (X + t(X))/2

  if (missing(phi)) #return starting values and boundaries for optimization
  {
    ones <- matrix(1, nrow(E), 1)
    Ed   <- diag(E)
    R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
    Rl   <- R[lower.tri(R)]
    Sl   <- S[lower.tri(S)]
    fit  <- nlme::gls(Sl ~ Rl, method = "ML")

    if (!nonnegative || coef(fit)[2] > 0) 
    {
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = phi[2], "tau" = -2 * log(sigma(fit)))
    }
    else
    {
      fit <- nlme::gls(Sl ~ 1, method = "ML")
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = 0, "tau" = -2 * log(sigma(fit)))
    }

    return(list(phi = phi, 
                lower = if(nonnegative) c(-Inf, 0, -Inf) else c(-Inf, -Inf, -Inf), 
                upper = c(Inf, Inf, Inf)))
  }
  else if (!(is.matrix(E)    & 
             is.matrix(S)    & 
             all(dim(E)  == dim(S)) &
             is.numeric(phi) & 
             length(phi) == 3 ))
    stop ("invalid inputs")

  names(phi) <- c("alpha", "beta", "tau")

  alpha <- phi["alpha"]
  tau   <- exp(phi["tau"])
  beta  <- phi["beta"]

  ones <- matrix(1, nrow(E), 1)
  Ed   <- diag(E)
  R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
  Rl   <- R[lower.tri(R)]
  Sl   <- S[lower.tri(S)]

  unos   <- matrix(1, length(Sl), 1)
  e      <- Sl - alpha * unos - beta * Rl
  loglik <- -0.5 * tau * t(e) %*% e + 0.5 * nrow(e) * log(tau)

  distance <- matrix(0, nrow(S), ncol(S))
  distance[lower.tri(distance)] <- Rl
  distance <- distance + t(distance)
  fitted <- alpha + beta * distance

  # gradients, hessians, mixed partial derivatives
  if (gradient || hessian || partial)
  {
    dPhi    <- matrix(0, length(phi), 1)
    ddPhi   <- matrix(0, length(phi), length(phi))
    ddEdPhi <- matrix(0, length(Rl),  length(phi))
    ddPhidS <- matrix(0, length(phi), length(Sl))
    rownames(dPhi) <- colnames(ddPhi) <- 
      rownames(ddPhi) <- colnames(ddEdPhi) <- 
        rownames(ddPhidS) <- names(phi)

    # gradient, phi
    dPhi["alpha",] <- t(unos) %*% e * tau
    dPhi["beta",]  <- t(Rl) %*% e * tau
    dPhi["tau",]   <- -0.5 * tau * t(e) %*% e + 0.5 * length(e)

    if (hessian || partial)
    {
      # hessian, phi x phi
      ddPhi["alpha", "alpha"] <- -t(unos) %*% unos * tau
      ddPhi["alpha",  "beta"] <- -tau * t(unos) %*% Rl
      ddPhi["alpha",   "tau"] <- tau * t(unos) %*% e
      ddPhi[ "beta",  "beta"] <- -t(Rl) %*% Rl * tau
      ddPhi[ "beta",   "tau"] <- tau * t(Rl) %*% e
      ddPhi[  "tau",   "tau"] <- -0.5 * tau * t(e) %*% e
      ddPhi                   <- ddPhi + t(ddPhi)
      diag(ddPhi)             <- diag(ddPhi)/2

      if (partial)
      {
        # gradient wrt E
        dR <- matrix(0, nrow(R), ncol(R))
        dR[lower.tri(dR)] <- 2 * beta * tau * as.vector(e)
        dR <- symm(dR)
        dE <- diag(nrow(R)) * (dR %*% ones %*% t(ones)) - dR

        # hessian offdiagonal, E x phi
        ddEdPhi[, "alpha"] <- -2 * beta * tau * unos
        ddEdPhi[,  "beta"] <- 2 * tau * (e - beta * Rl)
        ddEdPhi[,   "tau"] <- 2 * beta * tau * e
        ddEdPhi            <- apply(ddEdPhi, 2, function(x) { X <- matrix(0,nrow(E),ncol(E)); X[lower.tri(X)] <- x; X <- symm(X); diag(nrow(E)) * (X %*% ones %*% t(ones)) - X })

        # hessian offdiagonal, S x phi
        ddPhidS["alpha",] <- unos * tau
        ddPhidS["beta",]  <- Rl * tau
        ddPhidS["tau",]   <- -tau * e

        # jacobian products (label these properly)
        jacobian_E <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE <- -beta^2 * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          -ddEdE
        }

        jacobian_S <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE <- -beta * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          diag(ddEdE) <- 0
          -ddEdE
        }
      }
    }
  }

  if (validate)
  {
    num_gradient    <- numDeriv::grad(function(x) leastsquares(E = E, phi = x, S = S)$objective, phi)
    num_hessian     <- numDeriv::hessian(function(x) leastsquares(E = E, phi = x, S = S)$objective, phi)
    num_gradient_E  <- symm(matrix(numDeriv::grad(function(x) leastsquares(E = x, phi = phi, S = S)$objective, E), nrow(E), ncol(E)))
    num_partial_E   <- numDeriv::jacobian(function(x) leastsquares(E = E, phi = x, S = S)$gradient_E, phi)
    num_partial_S   <- numDeriv::jacobian(function(x) leastsquares(E = E, phi = phi, S = x)$gradient, S)[,lower.tri(S)]
    num_jacobian_E  <- function(X) matrix(c(X) %*% numDeriv::jacobian(function(x) leastsquares(E = x, phi = phi, S = S)$gradient_E, E), nrow(X), ncol(X))
    num_jacobian_S  <- function(X) matrix(c(X) %*% numDeriv::jacobian(function(x) leastsquares(E = E, phi = phi, S = x)$gradient_E, S), nrow(X), ncol(X))
  }

  list(objective        = -c(loglik), 
       fitted           = fitted,
       boundary         = nonnegative && beta == 0,
       gradient         = if(!gradient) NULL else -dPhi,
       hessian          = if(!hessian)  NULL else -ddPhi,
       gradient_E       = if(!partial)  NULL else -dE, 
       partial_E        = if(!partial)  NULL else -ddEdPhi,   # partial_E[i,k] is d(dl/dE_i)/dPhi_k where i is linearized matrix index
       partial_S        = if(!partial)  NULL else -ddPhidS,   # partial_S[i,k] is d(dl/dPhi_i)/dS_k where k is linearized matrix index
       jacobian_E       = if(!partial)  NULL else jacobian_E, # function mapping vectorized dg/dE to d(dg/dE)/dE
       jacobian_S       = if(!partial)  NULL else jacobian_S, # function mapping vectorized dg/dE to d(dg/dE)/dS
       num_gradient     = if(!validate) NULL else num_gradient,
       num_hessian      = if(!validate) NULL else num_hessian,
       num_gradient_E   = if(!validate) NULL else num_gradient_E,
       num_partial_E    = if(!validate) NULL else num_partial_E,
       num_partial_S    = if(!validate) NULL else num_partial_S,
       num_jacobian_E   = if(!validate) NULL else num_jacobian_E,
       num_jacobian_S   = if(!validate) NULL else num_jacobian_S)
}
class(leastsquares) <- "radish_measurement_model"

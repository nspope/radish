#' Generalized Wishart distance regression
#'
#' A function of class "measurement_model" that calculates likelihood,
#' gradient, hessian, and partial derivatives of nuisance parameters and the
#' Laplacian generalized inverse, using the generalized Wishart model of
#' McCullagh (2009) and Peterson et al (2019).
#'
#' @param E A submatrix of the generalized inverse of the graph Laplacian (alternatively, a covariance matrix)
#' @param S A matrix of observed genetic distances
#' @param phi Nuisance parameters (see details)
#' @param nu Number of genetic markers
#' @param gradient Compute gradient of negative loglikelihood wrt phi?
#' @param hessian Compute Hessian matrix of negative loglikelihood wrt phi?
#' @param partial Compute second partial derivatives of negative loglikelihood wrt phi and spatial covariates/observed genetic distances
#' @param nonnegative Ignored
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#'
#' @details The nuisance parameters are the scaling of the generalized inverse of the graph Laplacian (can be zero) and 
#' the log scalar multiple of the identity matrix that is added to the generalized inverse.
#'
#' TODO: formula
#'
#' @return A list containing at a minimum:
#'  \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
#' Additionally, if 'objective == TRUE':
#'  \item{objective}{(if 'objective') the negative loglikelihood}
#'  \item{fitted}{((if 'objective') a matrix of expected genetic distances among target vertices}
#'  \item{boundary}{(if 'objective') is the MLE on the boundary (e.g. no genetic structure)?}
#'  \item{gradient}{(if 'gradient') gradient of negative loglikelihood with respect to phi}
#'  \item{hessian}{(if 'hessian') Hessian matrix of the negative loglikelihood with respect to phi}
#'  \item{gradient_E}{(if 'partial') gradient with respect to the generalized inverse of the graph Laplacian}
#'  \item{partial_E}{(if 'partial') Jacobian of 'gradient_E' with respect to phi}
#'  \item{partial_S}{(if 'partial') Jacobian of 'gradient_E' with respect to phi}
#'  \item{jacobian_E}{(if 'partial') a function used for reverse algorithmic differentiation}
#'  \item{jacobian_S}{(if 'partial') a function used for reverse algorithmic differentiation}
#'
#' @references
#' McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica Sinica 19
#' Peterson et al. 2019. TODO
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#'
#' # null (IBD) resistance distance
#' laplacian_inv <- radish_distance(radish::loglinear_conductance, surface, 
#'                                  theta = matrix(0, 1, 2), covariance = TRUE)$covariance[,,1]
#' 
#' generalized_wishart(laplacian_inv, melip.Fst, phi = c(0.1, 0.5), nu = 1000)
#'
#' @export

generalized_wishart <- function(E, S, phi, nu, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  symm <- function(X) (X + t(X))/2

  if (missing(phi)) #return starting values and boundaries for optimization of phi
  {
    return(list(phi = c(1, 0), lower = c(0, -Inf), upper = c(Inf, Inf)))
  }
  else if (!(is.matrix(E)    & 
             is.matrix(S)    & 
             all(dim(E)  == dim(S)) &
             is.numeric(phi) & 
             length(phi) == 2 ))
    stop ("invalid inputs")

  stopifnot(nu > 0)

  # ensure that S is symmetric with 0 diagonal
  #TODO: warnings? check that all positive?
  S <- symm(S)
  diag(S) <- 0
  if (any(S < 0))
    warning("Some distances are negative")

  names(phi) <- c("tau", "sigma")
  tau   <- phi["tau"]
  sigma <- exp(phi["sigma"])

  # density is undefined if tau is negative
  stopifnot(tau >= 0)
  nonnegative <- TRUE

  ones      <- matrix(1, nrow(S), 1)
  I         <- diag(nrow(S))
  Sigma     <- tau * E + sigma * I
  SigInvOne <- solve(Sigma, ones)

  W        <- I - ones %*% solve(t(ones) %*% SigInvOne) %*% t(SigInvOne)
  SigInvW  <- solve(Sigma, W)
  eigSigW  <- eigen(SigInvW)
  P        <- eigSigW$vectors[,-nrow(Sigma)]
  D        <- diag(eigSigW$values[-nrow(Sigma)])
  ginvSigW <- P %*% solve(D) %*% t(P)

  loglik <- nu/4 * sum(diag(SigInvW %*% S)) + nu/2 * sum(log(diag(D)))

  fitted <- diag(Sigma) %*% t(ones) + ones %*% t(diag(Sigma)) - 2 * Sigma

  if (gradient || hessian || partial)
  {
    dPhi    <- matrix(0, length(phi), 1)
    ddPhi   <- matrix(0, length(phi), length(phi))
    ddEdPhi <- matrix(0, length(E),   length(phi))
    ddPhidS <- matrix(0, length(phi), sum(lower.tri(S)))
    rownames(dPhi) <- colnames(ddPhi) <- 
      rownames(ddPhi) <- colnames(ddEdPhi) <- 
        rownames(ddPhidS) <- names(phi)

    grad_Sigma <- -nu/2 * SigInvW - nu/4 * SigInvW %*% S %*% t(SigInvW)

    # gradient, phi
    dPhi["tau",]   <- sum(E * grad_Sigma) 
    dPhi["sigma",] <- sum(diag(grad_Sigma)) * sigma

    if (hessian || partial)
    {
      # see Golub GH, Pereyra V. 1973. The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate. SIAM Journal on Numerical Analysis 10(2): 413-432
      # ^^actually unnecessary
      dSigInvW_dtau <- -solve(Sigma, E %*% SigInvW)
      dSigInvW_dsigma <- -solve(Sigma, SigInvW)

      fuckoff_tau <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% E %*% solve(Sigma) +
        -c(t(ones) %*% solve(Sigma) %*% E %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
      #dgrad_dtau <- -0.5 * nu * dSigInvW_dtau - 0.25 * nu * dSigInvW_dtau %*% S %*% t(SigInvW) -
      #               0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dtau)
      dgrad_dtau <- -0.5 * nu * dSigInvW_dtau %*% W - 0.25 * nu * dSigInvW_dtau %*% S %*% t(SigInvW) -
                       0.5 * nu * W %*% dSigInvW_dtau - 0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dtau) +
                       0.5 * nu * W %*% dSigInvW_dtau %*% W
      dgrad_dtau <- -0.5 * nu * solve(Sigma) %*% fuckoff_tau - 
        nu/4 * solve(Sigma) %*% fuckoff_tau %*% S %*% t(W) %*% solve(Sigma) -
        nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_tau) %*% solve(Sigma) + dgrad_dtau
      
      fuckoff_sigma <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% solve(Sigma) +
        -c(t(ones) %*% solve(Sigma) %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
      dgrad_dsigma <- -0.5 * nu * dSigInvW_dsigma %*% W - 0.25 * nu * dSigInvW_dsigma %*% S %*% t(SigInvW) -
                       0.5 * nu * W %*% dSigInvW_dsigma - 0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dsigma) +
                       0.5 * nu * W %*% dSigInvW_dsigma %*% W
      dgrad_dsigma <- -0.5 * nu * solve(Sigma) %*% fuckoff_sigma - 
        nu/4 * solve(Sigma) %*% fuckoff_sigma %*% S %*% t(W) %*% solve(Sigma) -
        nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_sigma) %*% solve(Sigma) + dgrad_dsigma
      dgrad_dsigma <- dgrad_dsigma * sigma

      # hessian, phi x phi
      ddPhi["tau","tau"] <- sum(E * dgrad_dtau)
      ddPhi["tau","sigma"] <- sum(E * dgrad_dsigma)
      ddPhi["sigma","sigma"] <- sigma * sum(diag(dgrad_dsigma)) + dPhi["sigma",]
      ddPhi <- ddPhi + t(ddPhi)
      diag(ddPhi) <- diag(ddPhi)/2

      if(partial)
      {
        # gradient wrt E
        dE <- tau * grad_Sigma

        # hessian offdiagonal, E x phi
        ddEdPhi[,"tau"] <- grad_Sigma + tau * dgrad_dtau
        ddEdPhi[,"sigma"] <- dgrad_dsigma * tau

        # hessian offdiagonal, S x phi
        ddtaudS <- -nu/4 * SigInvW %*% E %*% t(SigInvW)
        ddsigmadS <- -nu/4 * SigInvW %*% t(SigInvW) * sigma
        ddPhidS["tau",] <- ddtaudS[lower.tri(ddtaudS)]
        ddPhidS["sigma",] <- ddsigmadS[lower.tri(ddsigmadS)]

        # jacobian products (label these properly)
        jacobian_E <- function(dE)
        {
          dE <- symm(dE)
          dSigInvW_dE <- -solve(Sigma, dE %*% SigInvW)
          fuckoff_E <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% dE %*% solve(Sigma) +
            -c(t(ones) %*% solve(Sigma) %*% dE %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
          dgrad_dE <- -0.5 * nu * dSigInvW_dE - 0.25 * nu * dSigInvW_dE %*% S %*% t(SigInvW) -
            0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dE)
          dgrad_dE <- -0.5 * nu * solve(Sigma) %*% fuckoff_E - 
            nu/4 * solve(Sigma) %*% fuckoff_E %*% S %*% t(W) %*% solve(Sigma) -
            nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_E) %*% solve(Sigma) + dgrad_dE
          -dgrad_dE * tau^2
        } #jesus christ

        jacobian_S <- function(dE)
        {
          dE <- symm(dE)
          dEdS <- -nu/4 * t(SigInvW) %*% dE %*% SigInvW * tau
          diag(dEdS) <- 0
          -dEdS
        } 
      }
    }
  }

  if (validate)
  {
    arg <- num_gradient <- numDeriv::jacobian(function(x) 
                                   generalized_wishart(E = E, 
                                        phi = x, 
                                        nu = nu,
                                        S = S)$grSig, 
                                   phi)
    num_gradient <- numDeriv::grad(function(x) 
                                   generalized_wishart(E = E, 
                                        phi = x, 
                                        nu = nu,
                                        S = S)$objective, 
                                   phi)

    num_hessian <- numDeriv::jacobian(function(x) 
                                     generalized_wishart(E = E, 
                                          phi = x, 
                                          nu = nu,
                                          S = S)$gradient, 
                                     phi)

    num_gradient_E <- symm(matrix(numDeriv::grad(function(x) 
                                                 generalized_wishart(E = x, 
                                                      phi = phi, 
                                                      nu = nu,
                                                      S = S)$objective, 
                                                 E), 
                                  nrow(E), ncol(E)))

    num_partial_E <- numDeriv::jacobian(function(x) 
                                        generalized_wishart(E = E, 
                                             phi = x, 
                                             nu = nu,
                                             S = S)$gradient_E, 
                                        phi)

    num_partial_S <- numDeriv::jacobian(function(x) 
                                        generalized_wishart(E = E, 
                                             phi = phi, 
                                             nu = nu,
                                             S = x)$gradient, 
                                        S)[,lower.tri(S)]

    num_jacobian_E <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         generalized_wishart(E = x, 
                                              phi = phi, 
                                              nu = nu,
                                              S = S)$gradient_E, 
                                         E), 
             nrow(X), ncol(X))

    num_jacobian_S <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         generalized_wishart(E = E, 
                                              phi = phi, 
                                              nu = nu,
                                              S = x)$gradient_E, 
                                         S), 
             nrow(X), ncol(X))
  }

  list(objective        = -c(loglik), 
       fitted           = fitted,
       boundary         = nonnegative && tau == 0,
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
class(generalized_wishart) <- "radish_measurement_model"

#' Maximum likelihood population effects
#'
#' A function of class "measurement_model" that calculates likelihood,
#' gradient, hessian, and partial derivatives of nuisance parameters and the
#' Laplacian generalized inverse, using the "maximum likelihood population
#' effects" model of Clarke et al (2002).
#'
#' @param E A submatrix of the generalized inverse of the graph Laplacian (alternatively, a covariance matrix)
#' @param S A matrix of observed genetic distances
#' @param nu Number of genetic markers (ignored)
#' @param phi Nuisance parameters (see details)
#' @param gradient Compute gradient of negative loglikelihood wrt phi?
#' @param hessian Compute Hessian matrix of negative loglikelihood wrt phi?
#' @param partial Compute second partial derivatives of negative loglikelihood wrt phi and spatial covariates/observed genetic distances
#' @param nonnegative Force slope to be nonnegative?
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#'
#' @details The nuisance parameters are the intercept, slope, log residual
#' standard deviation, and logit correlation parameter of the MLPE regression. If 'phi' is not
#' supplied, the MLE of phi is returned using package 'corMLPE' (github.com/nspope/corMLPE)
#'
#' TODO: formula
#'
#' @return A list containing at a minimum:
#'  \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
#' Additionally, if 'objective == TRUE':
#'  \item{objective}{the negative loglikelihood}
#'  \item{fitted}{matrix of expected genetic distances among target vertices}
#'  \item{boundary}{is the solution on the boundary (e.g. no genetic structure)?}
#'  \item{gradient}{(if 'gradient') gradient of negative loglikelihood with respect to phi}
#'  \item{hessian}{(if 'hessian') Hessian matrix of the negative loglikelihood with respect to phi}
#'  \item{gradient_E}{(if 'partial') gradient with respect to the generalized inverse of the graph Laplacian}
#'  \item{partial_E}{(if 'partial') Jacobian of 'gradient_E' with respect to phi}
#'  \item{partial_S}{(if 'partial') Jacobian of 'gradient_E' with respect to phi}
#'  \item{jacobian_E}{(if 'partial') a function used for reverse algorithmic differentiation}
#'  \item{jacobian_S}{(if 'partial') a function used for reverse algorithmic differentiation}
#'
#' @references
#' Clarke et al. TODO
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
#' mlpe(laplacian_inv, melip.Fst, nu = NULL) #without 'phi': return MLE of phi
#' mlpe(laplacian_inv, melip.Fst, nu = NULL, phi = c(0., 0.5, 0.1, 0.3))
#'
#' @export

mlpe <- function(E, S, nu, phi, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  symm <- function(X) (X + t(X))/2

  if (missing(phi)) #return starting values and boundaries for optimization of phi
  {
    ones <- matrix(1, nrow(E), 1)
    Ed   <- diag(E)
    R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
    Rl   <- R[lower.tri(R)]
    Sl   <- S[lower.tri(S)]
    Ind  <- which(lower.tri(R), arr.ind = TRUE)

    fit  <- nlme::gls(Sl ~ Rl, method = "ML", correlation = corMLPE::corMLPE(form = ~Ind1 + Ind2), 
                      data = data.frame(Sl, Rl, Ind1 = Ind[,1], Ind2 = Ind[,2]))

    if (!nonnegative || coef(fit)[2] > 0) 
    {
      rho <- fit$modelStruct$corStruct[[1]] #already unconstrained
      phi <- c("alpha" = coef(fit)[1], "beta" = coef(fit)[2], "tau" = -2 * log(sigma(fit)), 
               "rho" = rho)
    }
    else
    {
      fit  <- nlme::gls(Sl ~ 1, method = "ML", correlation = corMLPE::corMLPE(form = ~Ind1 + Ind2), 
                        data = data.frame(Sl, Ind1 = Ind[,1], Ind2 = Ind[,2]))
      rho <- fit$modelStruct$corStruct[[1]] #already unconstrained
      phi <- c("alpha" = coef(fit)[1], "beta" = 0, "tau" = -2 * log(sigma(fit)), 
               "rho" = rho)
    }

    return(list(phi   = phi, 
                lower = if (nonnegative) c(-Inf, 0, -Inf, -Inf) 
                        else c(-Inf, -Inf, -Inf, -Inf), 
                upper = c(Inf, Inf, Inf, Inf)))
  }
  else if (!(class(E)    == "matrix"    & 
             class(S)    == "matrix"    & 
             all(dim(E)  == dim(S))     &
             class(phi)  == "numeric"   & 
             length(phi) == 4           ))
    stop ("invalid inputs")

  names(phi) <- c("alpha", "beta", "tau", "rho")

  alpha <- phi["alpha"]
  beta  <- phi["beta"]
  tau   <- exp(phi["tau"])
  rho   <- 0.5 * plogis(phi["rho"])

  ones <- matrix(1, nrow(E), 1)
  Ed   <- diag(E)
  R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)

  Rl  <- R[lower.tri(R)]
  Sl  <- S[lower.tri(S)]
  Ind <- which(lower.tri(R), arr.ind = TRUE)

  unos   <- matrix(1, length(Sl), 1)
  U      <- Matrix::sparseMatrix(i = rep(1:length(Sl), 2), j = c(Ind), x = c(unos))

  if (nrow(E) > length(.mlpe_eigen))
  {
    warning("Cannot use precomputed eigendecomposition, consider updating `.mlpe_eigen`")
    eigUtU <- eigen(as.matrix(Matrix::t(U) %*% U)) #ideally this is precomputed
  }
  else
    eigUtU <- .mlpe_eigen[[nrow(E)]]
  D      <- eigUtU$values
  P      <- eigUtU$vectors
  Dr     <- D/(1-2*rho) + 1/rho

  SigmaInv <- function(x)
  {
    Ax <- 1/(1-2*rho) * x
    x <- Matrix::t(U) %*% Ax
    x <- t(P) %*% x
    x <- x / Dr
    x <- P %*% x
    x <- Ax - 1/(1-2*rho) * U %*% x
    as.matrix(x)
  }

  SigmaLogDet <- sum(log(Dr)) + length(D) * log(rho) + length(Sl) * log(1 - 2*rho)

  # products against inverse correlation matrix
  e       <- Sl - alpha * unos - beta * Rl
  Si_e    <- SigmaInv(e)
  loglik  <- -0.5 * tau * t(e) %*% Si_e + 0.5 * nrow(e) * log(tau) - 0.5 * SigmaLogDet 
  #check here that matches corMLPE

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

    drho_Si_e  <- Matrix::t(U) %*% Si_e
    drho_Si_e  <- as.matrix(2*Si_e - U %*% drho_Si_e)
    drho_trans <- rho * (1 - 2*rho)

    # gradient, phi
    dPhi["alpha",] <- t(unos) %*% Si_e * tau
    dPhi["beta",]  <- t(Rl) %*% Si_e * tau
    dPhi["tau",]   <- -0.5 * tau * t(e) %*% Si_e + 0.5 * length(e)
    dPhi["rho",]   <- 
      (-0.5 * tau * t(Si_e) %*% drho_Si_e - 
       0.5 * sum((2*D/(1-2*rho)^2 - 1/rho^2)/Dr) - 
       0.5 * length(D)/rho + length(Sl)/(1-2*rho)) * drho_trans

    if (hessian || partial)
    {
      Si_unos      <- SigmaInv(unos)
      Si_Rl        <- SigmaInv(Rl)
      Si_Sl        <- SigmaInv(Sl)
      Si_drho_Si_e <- SigmaInv(drho_Si_e)

      # hessian, phi x phi
      ddPhi["alpha", "alpha"] <- -tau * t(unos) %*% Si_unos
      ddPhi["alpha",  "beta"] <- -tau * t(unos) %*% Si_Rl
      ddPhi["alpha",   "tau"] <- tau * t(unos) %*% Si_e
      ddPhi["alpha",   "rho"] <- tau * t(Si_unos) %*% drho_Si_e * drho_trans
      ddPhi[ "beta",  "beta"] <- -tau * t(Rl) %*% Si_Rl
      ddPhi[ "beta",   "tau"] <- tau * t(Rl) %*% Si_e
      ddPhi[ "beta",   "rho"] <- tau * t(Si_Rl) %*% drho_Si_e * drho_trans
      ddPhi[  "tau",   "tau"] <- -0.5 * tau * t(e) %*% Si_e
      ddPhi[  "tau",   "rho"] <- -0.5 * tau * t(Si_e) %*% drho_Si_e * drho_trans
      ddPhi[  "rho",   "rho"] <- 
        (-tau * t(drho_Si_e) %*% Si_drho_Si_e +
         -0.5 * sum((8*D/(1-2*rho)^3 + 2/rho^3)/Dr) + 0.5 * sum((2*D/(1-2*rho)^2 - 1/rho^2)^2/Dr^2) + 
         0.5 * length(D)/rho^2 + 2 * length(Sl)/(1-2*rho)^2) * drho_trans^2 + dPhi["rho",] * (1 - 4*rho)
      ddPhi                   <- ddPhi + t(ddPhi)
      diag(ddPhi)             <- diag(ddPhi)/2

      if (partial)
      {
        # gradient wrt E
        dR <- matrix(0, nrow(R), ncol(R))
        dR[lower.tri(dR)] <- 2 * beta * tau * as.vector(Si_e)
        dR <- symm(dR)
        dE <- diag(nrow(R)) * (dR %*% ones %*% t(ones)) - dR

        # hessian offdiagonal, E x phi
        ddEdPhi[, "alpha"] <- -2 * beta * tau * Si_unos
        ddEdPhi[,  "beta"] <- 2 * tau * (Si_e - beta * Si_Rl)
        ddEdPhi[,   "tau"] <- 2 * beta * tau * Si_e
        ddEdPhi[,   "rho"] <- 2 * beta * tau * Si_drho_Si_e * drho_trans
        ddEdPhi            <- apply(ddEdPhi, 2, function(x) { X <- matrix(0,nrow(E),ncol(E)); X[lower.tri(X)] <- x; X <- symm(X); diag(nrow(E)) * (X %*% ones %*% t(ones)) - X })

        # hessian offdiagonal, S x phi
        ddPhidS["alpha",] <- tau * Si_unos
        ddPhidS["beta",]  <- tau * Si_Rl
        ddPhidS["tau",]   <- -tau * Si_e
        ddPhidS["rho",]   <- -tau * Si_drho_Si_e * drho_trans

        # jacobian products (label these properly)
        jacobian_E <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE[lower.tri(ddEdE)] <- SigmaInv(ddEdE[lower.tri(ddEdE)])
          ddEdE[upper.tri(ddEdE)] <- 0
          ddEdE <- ddEdE + t(ddEdE)
          ddEdE <- -beta^2 * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          -ddEdE
        }

        jacobian_S <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE[lower.tri(ddEdE)] <- SigmaInv(ddEdE[lower.tri(ddEdE)])
          ddEdE[upper.tri(ddEdE)] <- 0
          ddEdE <- ddEdE + t(ddEdE)
          ddEdE <- -beta * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          diag(ddEdE) <- 0
          -ddEdE
        }
      }
    }
  }

  # numerical validation
  if (validate)
  {
    num_gradient <- numDeriv::grad(function(x) 
                                   mlpe(E = E, 
                                        phi = x, 
                                        S = S)$objective, 
                                   phi)

    num_hessian <- numDeriv::hessian(function(x) 
                                     mlpe(E = E, 
                                          phi = x, 
                                          S = S)$objective, 
                                     phi)

    num_gradient_E <- symm(matrix(numDeriv::grad(function(x) 
                                                 mlpe(E = x, 
                                                      phi = phi, 
                                                      S = S)$objective, 
                                                 E), 
                                  nrow(E), ncol(E)))

    num_partial_E <- numDeriv::jacobian(function(x) 
                                        mlpe(E = E, 
                                             phi = x, 
                                             S = S)$gradient_E, 
                                        phi)

    num_partial_S <- numDeriv::jacobian(function(x) 
                                        mlpe(E = E, 
                                             phi = phi, 
                                             S = x)$gradient, 
                                        S)[,lower.tri(S)]

    num_jacobian_E <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         mlpe(E = x, 
                                              phi = phi, 
                                              S = S)$gradient_E, 
                                         E), 
             nrow(X), ncol(X))

    num_jacobian_S <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         mlpe(E = E, 
                                              phi = phi, 
                                              S = x)$gradient_E, 
                                         S), 
             nrow(X), ncol(X))
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
class(mlpe) <- "radish_measurement_model"

# pre-compute needed eigendecomposition for 3-100 demes
# eventually I could store this for larger range in /include?
.mlpe_eigen <- lapply(1:100, function(n) {
                        if (n == 1) NA
                        else {
                          Ind    <- which(lower.tri(diag(n)), arr.ind = TRUE)
                          U      <- Matrix::sparseMatrix(i = rep(1:nrow(Ind), 2), j = c(Ind), x = rep(1, nrow(Ind)*2))
                          eigen(as.matrix(Matrix::t(U) %*% U))
                        }
       })

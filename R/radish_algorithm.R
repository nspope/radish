#' Likelihood of parameterized conductance surface
#'
#' Calculates likelihood, gradient, hessian, and partial derivatives of a
#' parameterized conductance surface, given a function mapping spatial data to
#' conductance and a function mapping resistance distance (covariance) to
#' genetic distance; using the algorithm in Pope (in prep).
#'
#' @param f A function of class 'conductance_model'
#' @param g A function of class 'measurement_model'
#' @param s An object of class 'radish_graph'
#' @param S A matrix of observed genetic distances
#' @param theta Parameters for conductance surface (e.g. inputs to 'f')
#' @param nu Number of genetic markers (potentially used by 'g')
#' @param objective Compute negative loglikelihood?
#' @param gradient Compute gradient of negative loglikelihood wrt theta?
#' @param hessian Compute Hessian matrix of negative loglikelihood wrt theta?
#' @param partial Compute partial derivatives of negative loglikelihood wrt theta and spatial covariates/observed genetic distances
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#'
#' @return A list containing at a minimum:
#'  \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
#' Additionally, if 'objective == TRUE':
#'  \item{objective}{(if 'objective') the negative loglikelihood}
#'  \item{phi}{(if 'objective') fitted values of the nuisance parameters of 'g'}
#'  \item{boundary}{(if 'objective') is the solution on the boundary (e.g. no genetic structure)?}
#'  \item{fitted}{(if 'objective') matrix of expected genetic distances among target vertices}
#'  \item{gradient}{(if 'gradient') gradient of negative loglikelihood with respect to theta}
#'  \item{hessian}{(if 'hessian') Hessian matrix of the negative loglikelihood with respect to theta}
#'  \item{partial_X}{(if 'partial') Jacobian of the gradient with respect to the spatial covariates}
#'  \item{partial_S}{(if 'partial') Jacobian of the gradient with respect to the observed genetic distances}
#'
#' @references
#'
#' Pope NS. In prep. Fast gradient-based optimization of resistance surfaces.
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' radish_algorithm(radish::loglinear_conductance, radish::leastsquares, surface, 
#'                  ifelse(melip.Fst < 0, 0, melip.Fst), nu = 1000, theta = c(-0.3, 0.3))
#'
#' @export

radish_algorithm <- function(f, g, s, S, theta, nu = NULL, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  stopifnot(class(f)  ==  "radish_conductance_model")
  stopifnot(class(g)  ==  "radish_measurement_model")
  stopifnot(class(s)  ==  "radish_graph"            )

  stopifnot(is.matrix(S)     )
  stopifnot(is.numeric(theta))

  stopifnot(length(s$demes) == nrow(S)  )
  stopifnot(        nrow(S) == ncol(S)  )

  symm <- function(X) (X + t(X))/2

  # conductance
  C <- f(theta)

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- length(C$conductance)
  Q     <- s$laplacian
  Q@x[] <- -C$conductance[s$adj[1,]+1] - C$conductance[s$adj[2,]+1]

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  ones <- matrix(1, N, 1)
  v    <- sqrt(ones / N)
  Z    <- Matrix::Diagonal(N)[,s$demes]
  Zn   <- In %*% Z - (In %*% v) %*% (t(v) %*% Z)
  LQn  <- Matrix::update(s$choleski, Qn)
  G    <- Matrix::solve(LQn, Zn)
  tG   <- t(as.matrix(G))
  E    <- Matrix::t(Zn) %*% G 

  if (objective || gradient || hessian)
  {
    # measurement model
    subproblem <- radish_subproblem(g = g, E = as.matrix(E), S = S, nu = nu,
                                    nonnegative = nonnegative,
                                    control = NewtonRaphsonControl(verbose = FALSE, 
                                                                   ftol = 1e-10, 
                                                                   ctol = 1e-10))
    phi        <- subproblem$phi
    loglik     <- subproblem$loglikelihood

    # gradient calculation
    grad      <- rep(0, length(theta))
    hess      <- matrix(0, length(theta), length(theta))
    partial_X <- array(0, c(N, length(theta), length(theta)))
    partial_S <- array(0, c(nrow(S), ncol(S), length(theta)))
    if (gradient || hessian || partial)
    {
      dl_dE    <- subproblem$gradient 
      dl_dQnG  <- dl_dE %*% tG
      dl_dC    <- backpropagate_laplacian_to_conductance(dl_dQnG, tG, s$adj)
      for(k in 1:length(theta))
        grad[k] <- c(dl_dC) %*% C$df__dtheta(k)

      # hessian and mixed partial derivative calculations
      if (hessian || partial)
      {
        for (k in 1:length(theta))
        {
          dgrad__ddl_dC   <- C$df__dtheta(k) 
          dgrad__ddl_dQn  <- Matrix::forceSymmetric(backpropagate_conductance_to_laplacian(dgrad__ddl_dC, s$adj))
          dgrad__ddl_dQnG <- dgrad__ddl_dQn %*% G
          dgrad__ddl_dE   <- -tG %*% dgrad__ddl_dQnG
          dgrad__dE       <- subproblem$jacobian_E(as.matrix(dgrad__ddl_dE))
          dgrad__dG       <- Zn %*% dgrad__dE - 2 * dgrad__ddl_dQnG %*% dl_dE
          dgrad__dQnG     <- t(as.matrix(Matrix::solve(LQn, dgrad__dG)))
          dgrad__dC       <- backpropagate_laplacian_to_conductance(dgrad__dQnG, tG, s$adj)
          for(l in 1:length(theta))
            hess[k,l]       <- c(dgrad__dC) %*% C$df__dtheta(l) + c(dl_dC) %*% C$d2f__dtheta_dtheta(k, l)

          if (partial)
          {
            for(l in 1:length(theta))
              partial_X[,k,l] <- c(dgrad__dC) * C$df__dx(l) + c(dl_dC) * C$d2f__dtheta_dx(k, l)
            partial_S[,,k]  <- subproblem$jacobian_S(as.matrix(dgrad__ddl_dE))
          }
        }
      }
    }
  }

  # numerical validation
  if (validate)
  {
    #TODO the partial X won't work after formula refactor
    num_gradient <- numDeriv::grad(function(theta) 
                                   radish_algorithm(f = f,
                                                    g = g, 
                                                    s = s,
                                                    S = S, 
                                                    theta = theta)$objective, 
                                   theta, method = "Richardson")

    num_hessian  <- numDeriv::hessian(function(theta) 
                                      radish_algorithm(f = f,
                                                       g = g, 
                                                       s = s,
                                                       S = S, 
                                                       theta = theta)$objective, 
                                      theta, method = "Richardson")

    num_partial_X  <- array(0, c(N, length(theta), length(theta)))
    for (l in 1:length(theta))
      num_partial_X[,,l] <- t(numDeriv::jacobian(function(z) {
                                                 s$x[,l] <- z
                                                 radish_algorithm(f = f,
                                                                  g = g, 
                                                                  s = s,
                                                                  S = S, 
                                                                  theta = theta)$gradient },
                                                 s$x[,l], method="simple"))

    num_partial_S <- array(t(numDeriv::jacobian(function(S) 
                                                radish_algorithm(f = f,
                                                                 g = g, 
                                                                 s = s,
                                                                 S = S, 
                                                                 theta = theta)$gradient, 
                                                S, method = "simple")), 
                           c(nrow(S), ncol(S), length(theta)))
  }

  list (covariance    = E,
        objective     = if(!objective) NULL else loglik,
        phi           = if(!objective) NULL else phi,
        boundary      = if(!objective) NULL else subproblem$boundary, # the solution is on the boundary (e.g. no genetic structure) so all derivatives wrt theta are 0
        fitted        = if(!objective) NULL else subproblem$fit$fitted,
        gradient      = if(!gradient)  NULL else grad * (1 - subproblem$boundary), # wrt theta
        hessian       = if(!hessian)   NULL else hess * (1 - subproblem$boundary), # wrt theta
        partial_X     = if(!partial)   NULL else partial_X * (1 - subproblem$boundary), # partial_X[i,l,k] is \frac{\partial^2 L(theta,x)}{\partial theta_l \partial x_{ik}}
        partial_S     = if(!partial)   NULL else partial_S * (1 - subproblem$boundary), # partial_S[i,j,k] is \frac{\partial^2 L(theta,x)}{\partial theta_k \partial S_{ij}}
        num_gradient  = if(!validate)  NULL else num_gradient,
        num_hessian   = if(!validate)  NULL else num_hessian,
        num_partial_X = if(!validate)  NULL else num_partial_X,
        num_partial_S = if(!validate)  NULL else num_partial_S)
}

#' Likelihood of parameterized conductance surface
#'
#' Calculates likelihood, gradient, hessian, and partial derivatives
#' of the likelihood of a parameterized conductance surface,
#' given a function mapping spatial data to conductance and a function
#' mapping resistance distance (covariance) to genetic distance.
#'
#' @param f A function of class 'conductance_model'
#' @param g A function of class 'measurement_model'
#' @param s An object of class 'radish_graph'
#' @param S A matrix of observed genetic distances
#' @param theta Parameters for conductance surface (e.g. inputs to 'f')
#' @param objective Compute negative loglikelihood?
#' @param gradient Compute gradient of negative loglikelihood wrt theta?
#' @param hessian Compute Hessian matrix of negative loglikelihood wrt theta?
#' @param partial Compute partial derivatives of negative loglikelihood wrt theta and spatial covariates/observed genetic distances
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#'
#' @return A list containing
#'  \item{Something}{something}
#'  \item{Something}{something}
#'
#' @examples
#'  data(melip)
#'  make raster here
#'  surface <- radish_conductance_surface()
#'
#' @export

radish_algorithm <- function(f, g, s, S, theta = rep(0, ncol(s$x)), objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  stopifnot(    class(f) == "radish_conductance_model")
  stopifnot(    class(g) == "radish_measurement_model")
  stopifnot(    class(s) == "radish_graph"            )
  stopifnot(    class(S) == "matrix"                  )
  stopifnot(class(theta) == "numeric"                 )

  stopifnot(  length(theta) == ncol(s$x))
  stopifnot(length(s$demes) == nrow(S)  )
  stopifnot(        nrow(S) == ncol(S)  )

  symm <- function(X) (X + t(X))/2

  # conductance
  C <- f(s$x, theta)

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- nrow(s$x)
  Q     <-  s$laplacian
  Q@x[] <- -C$conductance[s$adj[1,]+1] - C$conductance[s$adj[2,]+1]

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
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
    subproblem <- radish_subproblem(g = g, E = as.matrix(E), S = S, 
                                    nonnegative = nonnegative,
                                    control = NewtonRaphsonControl(verbose = FALSE, 
                                                                   ftol = 1e-10, 
                                                                   ctol = 1e-10))
    phi        <- subproblem$phi
    loglik     <- subproblem$loglikelihood

    # gradient calculation
    grad      <- rep(0, length(theta))
    hess      <- matrix(0, length(theta), length(theta))
    partial_X <- array(0, c(nrow(s$x), length(theta), length(theta)))
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
    num_gradient <- numDeriv::grad(function(theta) 
                                   radish_algorithm(g = g, 
                                                    S = S, 
                                                    theta = theta, 
                                                    x = x, 
                                                    demes = demes, 
                                                    adj = adj)$objective, 
                                   theta, method = "Richardson")

    num_hessian  <- numDeriv::hessian(function(theta) 
                                      radish_algorithm(g = g, 
                                                       S = S, 
                                                       theta = theta, 
                                                       x = x, 
                                                       demes = demes, 
                                                       adj = adj)$objective, 
                                      theta, method = "Richardson")

    num_partial_X  <- array(0, c(nrow(s$x), length(theta), length(theta)))
    for (l in 1:length(theta))
      num_partial_X[,,l] <- t(numDeriv::jacobian(function(z) {
                                                 x[,l] <- z
                                                 radish_algorithm(g = g, 
                                                                  S = S, 
                                                                  theta = theta, 
                                                                  x = x, 
                                                                  demes = demes, 
                                                                  adj = adj)$gradient },
                                                  x[,l], method="simple"))

    num_partial_S <- array(t(numDeriv::jacobian(function(S) 
                                                radish_algorithm(g = g, 
                                                                 S = S, 
                                                                 theta = theta, 
                                                                 x = x, 
                                                                 demes = demes, 
                                                                 adj = adj)$gradient, 
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

# deprecated
DEPR_radish_algorithm <- function(f, g, s, S, theta = rep(0, ncol(s$x)), objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  stopifnot(    class(f) == "radish_conductance_model")
  stopifnot(    class(g) == "radish_measurement_model")
  stopifnot(    class(s) == "DEPR_radish_graph"       )
  stopifnot(    class(S) == "matrix"                  )
  stopifnot(class(theta) == "numeric"                 )

  stopifnot(  length(theta) == ncol(s$x))
  stopifnot(length(s$demes) == nrow(S)  )
  stopifnot(        nrow(S) == ncol(S)  )

  symm <- function(X) (X + t(X))/2

  # conductance
  C <- f(s$x, theta)

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N    <- nrow(s$x)
  Q    <- Matrix::sparseMatrix(i = s$adj[,1], j = s$adj[,2], dims = c(N, N),
                               x = -C$conductance[s$adj[,1]] - C$conductance[s$adj[,2]], 
                               use.last.ij = TRUE)
  Q    <- Matrix::forceSymmetric(Q)
  tadj <- rbind(Q@i, rep(1:Q@Dim[2] - 1, diff(Q@p))) #upper-triangular, 0-based

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  ones <- matrix(1, N, 1)
  v    <- sqrt(ones / N)
  Z    <- Matrix::Diagonal(N)[,s$demes]
  Zn   <- as.matrix(In %*% Z - (In %*% v) %*% (t(v) %*% Z))
  LQn  <- Matrix::Cholesky(Qn, LDL = TRUE)
  G    <- as.matrix(Matrix::solve(LQn, Zn)) #cnvt to dense
  E    <- t(Zn) %*% G #cnvt to dense

  if (objective || gradient || hessian)
  {
    # measurement model
    subproblem <- radish_subproblem(g = g, E = E, S = S, 
                                    nonnegative = nonnegative,
                                    control = NewtonRaphsonControl(verbose = FALSE, 
                                                                   ftol = 1e-10, 
                                                                   ctol = 1e-10))
    phi        <- subproblem$phi
    loglik     <- subproblem$loglikelihood

    # gradient calculation
    grad      <- rep(0, length(theta))
    hess      <- matrix(0, length(theta), length(theta))
    partial_X <- array(0, c(nrow(s$x), length(theta), length(theta)))
    partial_S <- array(0, c(nrow(S), ncol(S), length(theta)))
    if (gradient || hessian || partial)
    {
      dl_dE    <- subproblem$gradient 
      dl_dQnG  <- G %*% dl_dE
      dl_dC    <- backpropagate_laplacian_to_conductance(t(dl_dQnG), t(G), tadj)
      for(k in 1:length(theta))
        grad[k] <- c(dl_dC) %*% C$df__dtheta(k)

      # hessian and mixed partial derivative calculations
      if (hessian || partial)
      {
        for (k in 1:length(theta))
        {
          dgrad__ddl_dC   <- C$df__dtheta(k) 
          dgrad__ddl_dQn  <- Matrix::forceSymmetric(backpropagate_conductance_to_laplacian(dgrad__ddl_dC, tadj))
          dgrad__ddl_dQnG <- dgrad__ddl_dQn %*% G
          dgrad__ddl_dE   <- -as.matrix(t(G) %*% dgrad__ddl_dQnG) #cnvt to dense
          dgrad__dE       <- subproblem$jacobian_E(dgrad__ddl_dE)
          dgrad__dG       <- Zn %*% dgrad__dE - 2 * dgrad__ddl_dQnG %*% dl_dE
          dgrad__dQnG     <- as.matrix(Matrix::solve(LQn, dgrad__dG)) #cnvt to dense
          dgrad__dC       <- backpropagate_laplacian_to_conductance(t(dgrad__dQnG), t(G), tadj)
          for(l in 1:length(theta))
            hess[k,l]       <- c(dgrad__dC) %*% C$df__dtheta(l) + c(dl_dC) %*% C$d2f__dtheta_dtheta(k, l)

          if (partial)
          {
            for(l in 1:length(theta))
              partial_X[,k,l] <- c(dgrad__dC) * C$df__dx(l) + c(dl_dC) * C$d2f__dtheta_dx(k, l)
            partial_S[,,k]  <- subproblem$jacobian_S(dgrad__ddl_dE)
          }
        }
      }
    }
  }

  # numerical validation
  if (validate)
  {
    num_gradient <- numDeriv::grad(function(theta) 
                                   radish_algorithm(f = f, g = g, s = s, S = S, theta = theta)$objective, 
                                   theta, method = "Richardson")

    num_hessian  <- numDeriv::hessian(function(theta) 
                                      radish_algorithm(f = f, g = g, s = s, S = S, theta = theta)$objective, 
                                      theta, method = "Richardson")

    num_partial_X  <- array(0, c(nrow(s$x), length(theta), length(theta)))
    for (l in 1:length(theta))
      num_partial_X[,,l] <- t(numDeriv::jacobian(function(z) {
                                                   s$x[,l] <- z
                                                   radish_algorithm(f = f, g = g, s = s, S = S, theta = theta)$gradient },
                                                   s$x[,l], method="simple"))

    num_partial_S <- array(t(numDeriv::jacobian(function(S) 
                                                radish_algorithm(f = f, g = g, s = s, S = S, theta = theta)$gradient,
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

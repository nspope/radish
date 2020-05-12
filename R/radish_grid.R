#' Evaluate likelihood of a parameterized conductance surface
#'
#' Calculates profile likelihood of a parameterized conductance surface across
#' a grid of parameter values.
#'
#' @param f A function of class 'conductance_model'
#' @param g A function of class 'measurement_model'
#' @param s An object of class 'radish_graph'
#' @param S A matrix of observed genetic distances
#' @param nu Number of genetic markers (potentially used by 'g')
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param covariance If TRUE, additionally return (a submatrix of) the generalized inverse of graph Laplacian across the grid
#'
#' @return An object of class 'radish_grid'
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' theta <- as.matrix(expand.grid(x=seq(-6,6,length.out=21), y=seq(-6,6,length.out=21)))
#' grid <- radish_grid(radish::loglinear_conductance, radish::mlpe, surface, melip.Fst, theta, covariance=FALSE)
#' 
#' library(ggplot2)
#' ggplot(data.frame(loglik=grid$loglik, grid$theta)) + 
#'   geom_tile(aes(x=x,y=y,fill=-loglik)) + theme_bw() +
#'   geom_contour(aes(x=x,y=y,z=-loglik), color="black") +
#'   xlab(expression(theta[altitude])) +
#'   ylab(expression(theta[forestcover]))
#' 
#' @export

radish_grid <- function(f, g, s, S, theta, nu = NULL, nonnegative = TRUE, covariance = FALSE)
{
  stopifnot(is.matrix(theta))
  stopifnot(ncol(theta) == ncol(s$x))

  ll  <- rep(NA, nrow(theta))
  phi <- matrix(NA, length(g(S = S, E = rWishart(1, nrow(S), diag(nrow(S)))[,,1])$phi), nrow(theta))
  if (covariance)
    cv  <- array(NA, c(length(s$demes), length(s$demes), nrow(theta)))

  for (i in 1:nrow(theta))
    try({
      obj     <- radish_algorithm(f = f, g = g, s = s, S = S,
                                  theta = c(theta[i,]), nu = nu,
                                  gradient = FALSE,
                                  hessian  = FALSE,
                                  partial  = FALSE,
                                  nonnegative = nonnegative)
      ll[i]   <- obj$objective
      phi[,i] <- obj$phi
      if (covariance)
        cv[,,i] <- as.matrix(obj$covariance)
    })

  df <- list(theta = data.frame(theta), 
             loglik = ll,
             phi = phi,
             covariance = if(!covariance) NULL else cv)
  class(df) <- "radish_grid"
  df
}

#' Resistance distances from a parameterized conductance surface
#'
#' Calculates resistance distances associated with a parameterized conductance surface across
#' a grid of parameter values.
#'
#' @param f A function of class 'conductance_model'
#' @param s An object of class 'radish_graph'
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param covariance If TRUE, instead of resistance distances return the associated submatrix of the generalized inverse of graph Laplacian
#'
#' @return An object of class 'radish_grid'
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' theta <- as.matrix(expand.grid(x=seq(-6,6,length.out=21), y=seq(-6,6,length.out=21)))
#' distances <- radish_distance(radish::loglinear_conductance, surface, theta)
#' 
#' @export

radish_distance <- function(f, s, theta, covariance = FALSE)
{
  stopifnot(is.matrix(theta))
  stopifnot(ncol(theta) == ncol(s$x))

  cv <- array(NA, c(length(s$demes), length(s$demes), nrow(theta)))

  for (i in 1:nrow(theta))
    try({
      obj     <- radish_algorithm(f = f, g = leastsquares, s = s, S = diag(length(s$demes)), 
                                  theta = c(theta[i,]), 
                                  objective = FALSE,
                                  gradient  = FALSE,
                                  hessian   = FALSE,
                                  partial   = FALSE,
                                  nonnegative = nonnegative)
      cv[,,i] <- if(covariance) as.matrix(obj$covariance) else dist_from_cov(as.matrix(obj$covariance))
    })

  df <- list(theta = data.frame(theta))
  if (covariance) df$covariance <- cv else df$distance <- cv

  class(df) <- "radish_grid"
  df
}

#' Evaluate likelihood of a parameterized conductance surface
#'
#' Calculates profile likelihood of a parameterized conductance surface across
#' a grid of parameter values.
#'
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param formula A formula with the name of a matrix of observed genetic distances on the lhs, and covariates used in the creation of 'data' on the rhs
#' @param data An object of class 'radish_graph'
#' @param conductance_model A function of class 'radish_conductance_model_factory'
#' @param measurement_model A function of class 'radish_measurement_model'
#' @param nu Number of genetic markers (potentially used by 'g')
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param conductance If TRUE, edge conductance is the sum of cell conductances; otherwise edge conductance is the inverse of the sum of cell resistances (NOT USED; TODO)
#' @param covariance If TRUE, additionally return (a submatrix of) the generalized inverse of graph Laplacian across the grid
#'
#' @return An object of class 'radish_grid'
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=raster::scale(melip.altitude), forestcover=raster::scale(melip.forestcover)))
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' theta <- as.matrix(expand.grid(x=seq(-1,1,length.out=21), y=seq(-1,1,length.out=21)))
#' grid <- radish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
#'                     radish::loglinear_conductance, radish::mlpe)
#' 
#' library(ggplot2)
#' ggplot(data.frame(loglik=grid$loglik, grid$theta)) + 
#'   geom_tile(aes(x=x,y=y,fill=-loglik)) + theme_bw() +
#'   geom_contour(aes(x=x,y=y,z=-loglik), color="black") +
#'   xlab(expression(theta[altitude])) +
#'   ylab(expression(theta[forestcover]))
#' 
#' @export

radish_grid <- function(theta,
                        formula, 
                        data,
                        conductance_model = radish::loglinear_conductance, 
                        measurement_model = radish::mlpe, 
                        nu = NULL, 
                        nonnegative = TRUE, 
                        conductance = TRUE,
                        covariance  = FALSE)
{
  stopifnot(is.matrix(theta))

  # get response, remove lhs from formula
  terms    <- terms(formula)
  vars     <- as.character(attr(terms, "variables"))[-1]
  response <- attr(terms, "response")
  S        <- if(response) get(vars[attr(terms, "response")], parent.frame())
              else stop("'formula' must have genetic distance matrix on lhs")
  is_ibd   <- length(vars) == 1

  stopifnot(!is_ibd) #IBD; nothing to do

  formula  <- reformulate(attr(terms, "term.labels"))

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model <- conductance_model(formula, data$x) #TODO: is_ibd here and have factory modify accordingly
  default <- attr(conductance_model, "default")

  stopifnot(ncol(theta) == length(default))

  colnames(theta) <- names(default)

  ll  <- rep(NA, nrow(theta))
  phi <- matrix(NA, length(measurement_model(S = S, E = rWishart(1, nrow(S), diag(nrow(S)))[,,1])$phi), nrow(theta))
  if (covariance)
    cv  <- array(NA, c(length(data$demes), length(data$demes), nrow(theta)))

  for (i in 1:nrow(theta))
    try({
      obj     <- radish_algorithm(f = conductance_model, 
                                  g = measurement_model, 
                                  s = data, 
                                  S = S,
                                  theta = c(theta[i,]), nu = nu,
                                  gradient = FALSE,
                                  hessian  = FALSE,
                                  partial  = FALSE,
                                  nonnegative = nonnegative)
      ll[i]   <- -obj$objective
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
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param formula A formula with covariates used in the creation of 'data' on the rhs
#' @param data An object of class 'radish_graph'
#' @param conductance_model A function of class 'radish_conductance_model_factory'
#' @param conductance If TRUE, edge conductance is the sum of cell conductances; otherwise edge conductance is the inverse of the sum of cell resistances (NOT USED; TODO)
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

radish_distance <- function(theta,
                            formula, 
                            data,
                            conductance_model = radish::loglinear_conductance, 
                            conductance = TRUE,
                            covariance  = FALSE)
{
  stopifnot(is.matrix(theta))

  # get response, remove lhs from formula
  terms    <- terms(formula)
  is_ibd   <- length(attr(terms, "factors")) == 0

  stopifnot(!is_ibd) #IBD; nothing to do

  formula  <- reformulate(attr(terms, "term.labels"))

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model <- conductance_model(formula, data$x) #TODO: is_ibd here and have factory modify accordingly
  default <- attr(conductance_model, "default")

  stopifnot(ncol(theta) == length(default))

  colnames(theta) <- names(default)

  cv <- array(NA, c(length(data$demes), length(data$demes), nrow(theta)))

  for (i in 1:nrow(theta))
    try({
      obj     <- radish_algorithm(f = conductance_model, 
                                  g = leastsquares, #need something here, so...
                                  s = data, 
                                  S = diag(length(data$demes)), #need something here, so...
                                  theta = c(theta[i,]), 
                                  objective = FALSE,
                                  gradient  = FALSE,
                                  hessian   = FALSE,
                                  partial   = FALSE,
                                  nonnegative = TRUE)
      cv[,,i] <- if(covariance) as.matrix(obj$covariance) else dist_from_cov(as.matrix(obj$covariance))
    })

  df <- list(theta = data.frame(theta))
  if (covariance) df$covariance <- cv else df$distance <- cv

  class(df) <- "radish_grid"
  df
}

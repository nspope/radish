#' Log-link conductance model
#'
#' A function of class "conductance_model" that represents a log-linear
#' mapping from spatial covariates to conductance
#'
#' @details The model is of the form
#' 
#'   C_i = exp(x_{i1} * theta_{1} + x_{i2} * theta_{2} + ...)
#'
#' where "C_i" is the conductance of vertex "i", "x_ij" is the value of spatial
#' covariate "j" at vertex "i", and "theta_j" is the parameter associated with
#' covariate "j".
#'
#' @export

loglinear_conductance <- function(x, theta)
{
  conductance        <- as.vector(exp(x %*% theta))

  stopifnot(all(conductance > 0))

  # first- and second-order derivatives
  df__dx             <- function(k)    conductance * theta[k]
  df__dtheta         <- function(k)    conductance * x[,k]
  d2f__dtheta_dtheta <- function(k, l) conductance * x[,k] * x[,l]
  d2f__dtheta_dx     <- function(k, l) conductance * ((k==l) + x[,k] * theta[l])

  # asymptotic confidence intervals
  confint <- function(x, theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
  {
    scale <- match.arg(scale)
    cond_sd <- sqrt(rowSums((x %*% vcov) * x))
    ci <- log(conductance) + qnorm((1 - quantile)/2) * cond_sd %*% t(c(1, -1))
    colnames(ci) <- c("lower", "upper")
    attr(ci, "quantile") <- quantile 
    if (scale == "linpred") 
      return (ci)
    else if (scale == "conductance")
      return (exp(ci))
  }

  list(conductance        = conductance,
       confint            = confint,
       df__dx             = df__dx,
       df__dtheta         = df__dtheta,
       d2f__dtheta_dtheta = d2f__dtheta_dtheta, 
       d2f__dtheta_dx     = d2f__dtheta_dx)
}
class(loglinear_conductance) <- c("radish_conductance_model")

#' Identity-link conductance model
#'
#' A function of class "conductance_model" that represents a linear
#' mapping from spatial covariates to conductance
#'
#' @details The model is of the form:
#'
#'   C_i = x_{i1} * theta_{1} + x_{i2} * theta_{2} + ...
#'
#' where "C_i" is the conductance of vertex "i", "x_ij" is the value of spatial
#' covariate "j" at vertex "i", and "theta_j" is the parameter associated with
#' covariate "j".
#'
#' @export

linear_conductance <- function(x, theta)
{
  conductance        <- as.vector(x %*% theta)

  stopifnot(all(conductance > 0))

  ones <- matrix(1, nrow(x), 1)

  # asymptotic confidence intervals
  confint <- function(x, theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
  {
    scale <- match.arg(scale)
    cond_sd <- sqrt(rowSums((x %*% vcov) * x))
    ci <- conductance + qnorm((1 - quantile)/2) * cond_sd %*% t(c(1, -1))
    colnames(ci) <- c("lower", "upper")
    attr(ci, "quantile") <- quantile 
    if (scale == "linpred") 
      return (ci)
    else if (scale == "conductance")
      return (ci)
  }

  # first- and second-order derivatives
  df__dx             <- function(k)    ones * theta[k]
  df__dtheta         <- function(k)    x[,k]
  d2f__dtheta_dtheta <- function(k, l) 0. * ones
  d2f__dtheta_dx     <- function(k, l) (k==l) * ones

  list(conductance        = conductance,
       confint            = confint,
       df__dx             = df__dx,
       df__dtheta         = df__dtheta,
       d2f__dtheta_dtheta = d2f__dtheta_dtheta, 
       d2f__dtheta_dx     = d2f__dtheta_dx)
}
class(linear_conductance) <- c("radish_conductance_model")

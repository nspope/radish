assemble_model_matrix <- function(formula, spdat)
{
  stopifnot(class(formula) == "formula")
  stopifnot(class(spdat) == "data.frame")

  # check if formula is consistant with data, remove response, add intercept
  formula_covariates <- attr(delete.response(terms(formula)), "factors")
  if (length(formula_covariates) > 0)
  {
    stopifnot(rownames(formula_covariates) %in% colnames(spdat))
    formula <- reformulate(colnames(formula_covariates))

    # if any layers are not in formula, remove them
    missing_covariates <- !(colnames(spdat) %in% rownames(formula_covariates))
    if (any(missing_covariates))
    {
      unused_covariates <- colnames(spdat)[missing_covariates]
      warning("Removed unused spatial covariates: ", 
              paste(unused_covariates, collapse = " "))
      spdat <- spdat[,!missing_covariates,drop=FALSE]
    }
  }
  else
    formula <- formula(~1)

  # get model matrix and check for rank deficiency
  # NOTE: sparse via Matrix::sparse.model.matrix?
  spdat <- model.matrix(formula, data = spdat)
  stopifnot(qr(spdat)$rank == ncol(spdat))
  if (ncol(spdat) > 1) #unless IBD, remove intercept
    spdat <- spdat[,colnames(spdat) != "(Intercept)", drop=FALSE]

  spdat
}

#' Log-link conductance model
#'
#' Returns a function of class "conductance_model" that represents a log-linear
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

loglinear_conductance <- function(formula, x)
{
  x <- assemble_model_matrix(formula, x)

  # default starting values
  default <- rep(0, ncol(x))
  names(default) <- colnames(x)

  conductance_model <- function(theta)
  {
    stopifnot(length(theta) == ncol(x))

    conductance        <- as.vector(exp(x %*% theta))

    stopifnot(all(conductance > 0))

    # first- and second-order derivatives
    df__dx             <- function(k)    conductance * theta[k]
    df__dtheta         <- function(k)    conductance * x[,k]
    d2f__dtheta_dtheta <- function(k, l) conductance * x[,k] * x[,l]
    d2f__dtheta_dx     <- function(k, l) conductance * ((k==l) + x[,k] * theta[l])

    # asymptotic confidence intervals
    confint <- function(theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
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

    # default starting values
    default <- function()
    {
      out <- rep(0, ncol(x))
      names(out) <- colnames(x)
      out
    }

    list(conductance        = conductance,
         confint            = confint,
         df__dx             = df__dx,
         df__dtheta         = df__dtheta,
         d2f__dtheta_dtheta = d2f__dtheta_dtheta, 
         d2f__dtheta_dx     = d2f__dtheta_dx)
  }

  class(conductance_model) <- c("radish_conductance_model")
  attr(conductance_model, "default") <- default
  conductance_model
}
class(loglinear_conductance) <- c("radish_conductance_model_factory")

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

linear_conductance <- function(formula, x)
{
  x <- assemble_model_matrix(formula, x)

  # default starting values
  default <- rep(1, ncol(x))
  names(default) <- colnames(x)

  conductance_model <- function(theta)
  {
    stopifnot(length(theta) == ncol(x))

    conductance        <- as.vector(x %*% theta)

    stopifnot(all(conductance > 0))

    ones <- matrix(1, nrow(x), 1)

    # asymptotic confidence intervals
    confint <- function(theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
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

  class(conductance_model) <- c("radish_conductance_model")
  attr(conductance_model, "default") <- default
  conductance_model
}
class(linear_conductance) <- c("radish_conductance_model_factory")

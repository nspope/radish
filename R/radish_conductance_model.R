loglinear_conductance <- function(x, theta)
{
  conductance        <- as.vector(exp(x %*% theta))

  stopifnot(all(conductance > 0))

  # first- and second-order derivatives
  df__dx             <- function(k)    conductance * theta[k]
  df__dtheta         <- function(k)    conductance * x[,k]
  d2f__dtheta_dtheta <- function(k, l) conductance * x[,k] * x[,l]
  d2f__dtheta_dx     <- function(k, l) conductance * ((k==l) + x[,k] * theta[l])

  list(conductance        = conductance,
       df__dx             = df__dx,
       df__dtheta         = df__dtheta,
       d2f__dtheta_dtheta = d2f__dtheta_dtheta, 
       d2f__dtheta_dx     = d2f__dtheta_dx)
}
class(loglinear_conductance) <- c("radish_conductance_model")

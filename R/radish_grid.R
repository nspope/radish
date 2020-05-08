radish_grid <- function(f, g, s, S, theta, nonnegative = TRUE, covariance = FALSE)
{
  stopifnot(class(theta) == "matrix" ) 
  stopifnot( ncol(theta) == ncol(s$x))

  ll  <- rep(NA, nrow(theta))
  phi <- matrix(NA, length(g(S = S, E = rWishart(1, nrow(S), diag(nrow(S)))[,,1])), nrow(theta))
  if (covariance)
    cv  <- array(NA, c(length(s$demes), length(s$demes), nrow(theta)))

  for (i in 1:nrow(theta))
    try({
      obj     <- radish_algorithm(f = f, g = g, s = s, S = S, 
                                  theta = c(theta[i,]), 
                                  gradient = FALSE,
                                  hessian  = FALSE,
                                  partial  = FALSE,
                                  nonnegative = nonnegative)
      ll[i]   <- obj$objective
      phi[,i] <- obj$phi
      if (covariance)
        cv[,,i] <- obj$covariance
    })

  df <- list(theta = data.frame(theta), 
             loglik = ll,
             phi = phi,
             covariance = if(!covariance) NULL else cv)
  class(df) <- "radish_grid"
  df
}

radish_distance <- function(f, s, theta, covariance = FALSE)
{
  stopifnot(class(theta) == "matrix" ) 
  stopifnot( ncol(theta) == ncol(s$x))

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
      cv[,,i] <- if(covariance) obj$covariance else dist_from_cov(obj$covariance)
    })

  df <- list(theta = data.frame(theta))
  if (covariance) df$covariance <- cv else df$distance <- cv

  class(df) <- "radish_grid"
  df
}

# simple backtracking to find finite function values
Backtracking <- function (dphifn, phi_0, dphi_0)
{
  maxit <- 10
  rho <- 1e-4
  c2 <- 0.9
  delta <- 0.2
  alpha <- 5
  for (iter in 1:maxit)
  {
    ev <- dphifn(alpha)
    val <- ev$objective
    gra <- ev$gradient
    if (is.finite(val) && is.finite(gra) && val <= phi_0)# && val <= phi_0 + alpha*rho*dphi_0 && -gra <= -c2*dphi_0)
    {
      return(alpha)
    }
    else
    {
      alpha <- alpha * delta
    }
  }
  warning("Maxit reached in backtracking")
  return(alpha)
}

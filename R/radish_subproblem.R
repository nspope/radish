radish_subproblem <- function(g, E, S, nu, phi = g(E = E, S = S, nonnegative = nonnegative), nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10, verbose = TRUE))
{
  # use Newton-Raphson to profile out nuisance parameters
  subproblem  <- BoxConstrainedNewton(phi$phi, 
                               function(par, gradient, hessian) 
                                 g(E = E, S = S, nu = nu, phi = c(par), 
                                   gradient = gradient, 
                                   hessian = hessian, 
                                   partial = FALSE, 
                                   nonnegative = nonnegative),
                               lower = phi$lower,
                               upper = phi$upper,
                               control = control)

  # refit, computing partial derivatives
  phi         <- subproblem$par
  fit         <- g(E = E, S = S, nu = nu, phi = c(phi), partial = TRUE, nonnegative = nonnegative)
  gradient_E  <- fit$gradient_E

  # for hessian, need to get d(dg/dE)/dE via adjoint method,
  #   dg/dE = \partial (dg/dE)/\partial E + \partial (dg/dE)/\partial \hat{phi} \times \partial \hat{\phi}/\partial E
  # where
  #   dg/dphi = 0 ==> d(dg/dphi)/dE = 0 ==> 
  #       \partial (dg/dphi)/\partial E + \partial (dg/dphi)/\partial phi \times dphi/dE = 0 ==>
  #       dphi/dE = -[\partial (dg/dphi)/\partial phi]^-1 \partial (dg/dphi)/\partial E
  partial_E   <- fit$partial_E
  invhess     <- MASS::ginv(fit$hessian)
  jacobian_E  <- function(dotdotE)
  { 
    #why is this nonzero when on boundary?
    dphi_dE   <- -matrix(c(dotdotE) %*% partial_E %*% invhess %*% t(partial_E), nrow(dotdotE), ncol(dotdotE))
    return (fit$jacobian_E(dotdotE) + dphi_dE)
  }

  # likewise, to get the leverage dtheta/dy, 
  #    dl/dtheta = 0 ==> d(dl/dtheta)/dy = 0 ==>
  #       \partial (dl/dtheta)/\partial y + \partial (dl/dtheta)/\partial theta \times dtheta/dy = 0 ==>
  #       dtheta/dy = -[\partial (dl/dtheta)/\partial theta]^{-1} \partial (dl/dtheta)/partial y
  # so, need the change in the gradient with y
  partial_S   <- fit$partial_S
  jacobian_S  <- function(dotdotE)
  { 
    #why is this nonzero when on boundary?
    dphi_dS                     <- matrix(0, nrow(S), ncol(S))
    dphi_dS[lower.tri(dphi_dS)] <- -c(dotdotE) %*% partial_E %*% invhess %*% partial_S
    dphi_dS                     <- dphi_dS + t(dphi_dS)
    return (fit$jacobian_S(dotdotE) + dphi_dS)
  }

  # numerical validation
  if (validate)
  {
    silence <- function(control) { control$verbose = FALSE; control }

    num_jacobian_E <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         radish_subproblem(g = g, 
                                                           E = x, 
                                                           S = S, 
                                                           phi = phi, 
                                                           nonnegative = nonnegative,
                                                           control = silence(control))$gradient, 
                                         E), 
             nrow(E), ncol(E))

    num_jacobian_S <- function(X) 
      matrix(c(X) %*% numDeriv::jacobian(function(x) 
                                         radish_subproblem(g = g, 
                                                           E = E, 
                                                           S = x, 
                                                           phi = phi, 
                                                           control = silence(control))$gradient, 
                                         S), 
             nrow(S), ncol(S))
  }

  list(fit            = fit,
       loglikelihood  = fit$objective,
       boundary       = fit$boundary,
       phi            = phi,
       gradient       = gradient_E,
       jacobian_E     = jacobian_E,
       jacobian_S     = jacobian_S,
       num_jacobian_E = if(!validate) NULL else num_jacobian_E,
       num_jacobian_S = if(!validate) NULL else num_jacobian_S)
}

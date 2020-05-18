#' Controls for Newton-like optimizer
#'
#' Control settings for Newton and quasi-Newton algorithms
#'
#' @param maxit Maximum number of Newton steps
#' @param ctol Convergence tolerance in the gradient (max(abs(grad)) < ctol)
#' @param ftol Convergence tolerance in the objective (f - fold < ftol)
#' @param etol Maximum spectral gap for a eigenvalue to be considered too close to singular
#' @param verbose Print progress to stdout
#' @param eps For box-constrained optimization, how close to a boundary can be considered "on the boundary"
#' @param del For quasi-Newton methods, a typical step size used to initialize the approximate Hessian
#' @param ls.control Control settings for the line search
#'
#' @export
NewtonRaphsonControl <- function(maxit = 100, 
                                 ctol = sqrt(.Machine$double.eps), 
                                 ftol = sqrt(.Machine$double.eps), 
                                 etol = 10*.Machine$double.eps, 
                                 verbose = FALSE, 
                                 eps = 1e-8, 
                                 del = 1,
                                 ls.control = HagerZhangControl())
  list(maxit = maxit, ctol = ctol, etol = etol, 
       ftol = ftol, verbose = verbose, eps = eps, 
       del = del, ls.control = ls.control)

BoxConstrainedNewton <- function(par, fn, lower = rep(-Inf, length(par)), upper = rep(Inf, length(par)), control = NewtonRaphsonControl())
{
  BoxConstrainedNewtonNaN <- function()
  {
    list(objective = NaN,
         gradient  = matrix(NaN, length(par), 1),
         hessian   = matrix(NaN, length(par), length(par)))
  }
  
  prettify <- function(x)
    formatC(x, digits=3, width=5, format="e")

  zero_bounded_variables <- function(gradient, par, lower, upper, eps = 1e-8)
  {
    # set gradient to 0 for active constraints
    tol <- eps * abs(par)
    gradient <- ifelse(upper - tol <= par & gradient < 0, 0, gradient)
    gradient <- ifelse(lower + tol >= par & gradient > 0, 0, gradient)
    gradient
  }

  gap_step_bounded_variables <- function(desc, par, gradient, lower, upper, eps = 1e-8)
  {
    # modify search direction so that at alpha == 1, actively constrained variables are set to the boundary
    tol <- eps * abs(par)
    desc <- ifelse(upper - tol <= par & gradient < 0, upper - par, desc)
    desc <- ifelse(lower + tol >= par & gradient > 0, lower - par, desc)
    desc
  }

  project <- function(x, lower, upper)
    pmin(pmax(x, lower), upper)

  stopifnot(lower < upper)

  for (var in names(control)) assign(var, control[[var]])
  etol <- etol * length(par)

  if (verbose)
    cat("Projected Newton-Raphson with Hager-Zhang line search\n")

  convergence <- 0
  par <- as.matrix(par)

  for (i in 1:maxit)
  {
    fit   <- fn(par, gradient = TRUE, hessian = TRUE)
    delta <- if (i > 1) abs(oldfit$objective - fit$objective) else 0

    if (verbose)
      cat(paste0("[", i, "]"), 
          "f(x) =", prettify(-fit$objective),
          "  |f(x)-fold(x)| =", prettify(delta),
          "  max|f'(x)| =", prettify(max(abs(fit$gradient))),
          "  |f''(x)| =", prettify(-det(fit$hessian)),
          "\n")

    if (max(abs(fit$gradient)) < ctol || (i > 1 && delta < ftol))
      break

    gradient     <- fit$gradient
    gradient_box <- zero_bounded_variables(gradient, par, lower, upper, eps)
    ehess        <- eigen(fit$hessian)
    ehess$values <- abs(ehess$values)
    ehess$values <- ifelse(ehess$values < max(abs(fit$hessian)) * etol, 1, ehess$values)
    ihess        <- ehess$vectors %*% solve(diag(ehess$values, nrow=length(par))) %*% t(ehess$vectors)
    desc         <- gap_step_bounded_variables(-ihess %*% gradient_box, par, gradient, lower, upper, eps)
    phi0         <- fit$objective
    dphi0        <- c(t(desc) %*% gradient_box)

    dphi_fn <- function(alpha) 
    {
      tryCatch({
        phi <- fn(project(par + alpha*desc, lower, upper), gradient = TRUE, hessian = FALSE)
        grb <- zero_bounded_variables(phi$gradient, par + alpha*desc, lower, upper, eps)
        list(objective = phi$objective, gradient = c(t(desc) %*% grb))
      }, error = function(e) {
        BoxConstrainedNewtonNaN()
      })
    }

    #alpha <- HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    alpha <- tryCatch({
      HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    }, error = function(err) {
      cat("... switched to backtracking\n")
      Backtracking(dphi_fn, phi0, dphi0)
    })
    par <- project(par + alpha*desc, lower, upper)

    oldfit <- fit
  }

  boundary_fit <- any(par == lower | par == upper)
  if (verbose)
    if (boundary_fit)
      cat ("Solution on boundary with `max(abs(gradient))` ==", max(abs(fit$gradient)), "and `diff(f)` ==", delta, "\n")
    else
      cat ("Solution on interior with `max(abs(gradient))` ==", max(abs(fit$gradient)), "and `diff(f)` ==", delta, "\n")

  if (i == maxit)
  {
    warning("`maxit` reached for Newton steps")
    convergence = 1
  } 

  list(par = par,
       gradient = fit$gradient,
       hessian = fit$hessian,
       value = fit$objective,
       fit = fit,
       iters = i,
       boundary = boundary_fit,
       convergence = convergence)
}

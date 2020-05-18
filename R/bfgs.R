#for now focus on stable BFGS
#setRefClass("LBFGSstorage", fields = list(ss = "matrix", yy = "matrix", m = "integer", k = "integer"))

BoxConstrainedBFGS <- function(par, fn, lower = rep(-Inf, length(par)), upper = rep(Inf, length(par)), control = NewtonRaphsonControl())
{
  BFGSNaN <- function()
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
    cat("BFGS with Hager-Zhang line search\n")

  convergence <- 0
  initialized <- 0
  par <- as.matrix(par)

  for (i in 1:maxit)
  {
    fit   <- fn(par, gradient = TRUE, hessian = FALSE)
    delta <- if (i > 1) abs(oldfit$objective - fit$objective) else 0

    if (verbose)
      cat(paste0("[", i, "]"), 
          "f(x) =", prettify(-fit$objective),
          "  |f(x) - fold(x)| =", prettify(delta),
          "  max|f'(x)| =", prettify(max(abs(fit$gradient))),
          "\n")

    if (max(abs(fit$gradient)) < ctol || (i > 1 && delta < ftol))
      break

    gradient     <- fit$gradient
    gradient_box <- zero_bounded_variables(gradient, par, lower, upper, eps)

    if(initialized > 0) 
    { #BFGS update from Nodecal and Wright Ch 6
      #with damping from Nodecal and Wright Ch 18 to ensure matrix is sufficiently positive definite
      yy    <- gradient - oldfit$gradient
      ss    <- alpha*desc
      irho  <- c(t(yy) %*% ss) 
      if(initialized == 1)
      { #rescale initial Hessian as per Nodecal and Wright 6.20
        initialized <- 2
        ihess       <- diag(nrow(ihess)) * irho/c(t(yy) %*% yy)
        hess        <- diag(nrow(hess)) * c(t(yy) %*% yy)/irho
      } 
      sBs   <- c(t(ss) %*% hess %*% ss)
      theta <- if (irho >= 0.2 * sBs) 1.0 else (0.8 * sBs)/(sBs - irho)
      if (verbose && theta < 1.0)
        cat("... damped BFGS update\n")
      rr    <- theta * yy + (1 - theta) * hess %*% ss
      rho   <- 1./c(t(rr) %*% ss) 
      upd   <- diag(nrow(ihess)) - rho * ss %*% t(rr)
      ihess <- upd %*% ihess %*% t(upd) + rho * ss %*% t(ss)
      hess  <- hess - hess %*% ss %*% t(ss) %*% hess / sBs + rho * rr %*% t(rr)
    } 
    else
    { #initial (diagonal) Hessian approximation
      initialized <- 1 
      ihess       <- del/sqrt(sum(gradient * gradient)) * diag(length(par))
      hess        <- sqrt(sum(gradient * gradient))/del * diag(length(par))
    }

    desc  <- gap_step_bounded_variables(-ihess %*% gradient_box, par, gradient, lower, upper, eps)
    phi0  <- fit$objective
    dphi0 <- c(t(desc) %*% gradient_box)

    dphi_fn <- function(alpha) 
    {
      tryCatch({
        phi <- fn(project(par + alpha*desc, lower, upper), gradient = TRUE, hessian = FALSE)
        grb <- zero_bounded_variables(phi$gradient, par + alpha*desc, lower, upper, eps)
        list(objective = phi$objective, gradient = c(t(desc) %*% grb))
      }, error = function(e) {
        # TODO: could "restart" hessian approximation here?
        # Not sure if this is necessary with damping
        BFGSNaN()
      })
    }

    alpha <- tryCatch({
      HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    }, error = function(err) {
      cat("Switched to backtracking\n")
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
    warning("`maxit` reached for quasi-Newton steps")
    convergence = 1
  } 

  list(par = par,
       gradient = fit$gradient,
       hessian = fit$hessian,
       value = fit$objective,
       fit = fit,
       boundary = boundary_fit,
       convergence = convergence)
}


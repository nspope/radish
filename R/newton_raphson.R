BFGS <- function(par, fn, lower = rep(-Inf, length(par)), upper = rep(Inf, length(par)), control = NewtonRaphsonControl())
{
  #TODO: this isn't finished

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
  par <- matrix(par, length(par), 1)

  for (i in 1:maxit)
  {
    fit   <- fn(par, gradient = TRUE, hessian = FALSE)
    delta <- if (i > 1) abs(oldfit$objective - fit$objective) else 0

    if (verbose)
      cat(paste0("[", i, "]"), 
          "f(x) =", prettify(-fit$objective),
          "  |\u0394 f(x)| =", prettify(delta),
          "  max|\u2207 f(x)| =", prettify(max(abs(fit$gradient))),
          "\n")

    if (max(abs(fit$gradient)) < ctol || (i > 1 && delta < ftol))
      break

    gradient     <- fit$gradient
    if(i > 1)
    {
      oldgrad <- oldfit$gradient
      yy <- gradient - oldgrad
      ss <- alpha*desc
      hess <- hess + yy %*% t(yy)/c(t(yy) %*% ss) - hess %*% ss %*% t(ss) %*% t(hess) / c(t(ss) %*% hess %*% ss)
    } else
    {
      hess <- sqrt(gradient * gradient) * diag(length(par))
    }

    desc         <- -solve(hess, gradient)
    phi0         <- fit$objective
    dphi0        <- c(t(desc) %*% gradient)

    dphi_fn <- function(alpha) 
    {
      tryCatch({
        phi <- fn(par + alpha*desc, gradient = TRUE, hessian = FALSE)
        grb <- phi$gradient
        list(objective = phi$objective, gradient = c(t(desc) %*% grb))
      }, error = function(e) {
        BFGSNaN()
      })
    }

    alpha <- tryCatch({
      HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    }, error = function(err) {
      cat("Switched to backtracking\n")
      Backtracking(dphi_fn, phi0, dphi0)
    })
    par   <- par + alpha*desc

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
       boundary = boundary_fit,
       convergence = convergence)

  # I don't understand why Armijo linesearch is needed
}

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
  par <- matrix(par, length(par), 1)

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
      cat("Switched to backtracking for remainder of step\n")
    Backtracking(dphi_fn, phi0, dphi0)
    })
    par   <- project(par + alpha*desc, lower, upper)

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

  # I don't understand why Armijo linesearch is needed
}

#' Control settings for Hager-Zhang line search
#'
#' TODO after finalizing line search
#'
#' @export
HagerZhangControl <- function(delta = 0.1, sigma = 0.9, alphamax = Inf, rho = 5.0, epsilon = 1e-6, gamma = 0.66, linesearchmax = 50, psi3 = 0.1, c = 1.0, verbose = FALSE)
  list(delta=delta, sigma=sigma, alphamax=alphamax, rho=rho, epsilon=epsilon, gamma=gamma, linesearchmax = linesearchmax, psi3 = psi3, c = c, verbose = verbose)

setRefClass("HagerZhangStorage", fields=list(alphas="numeric", values="numeric", slopes="numeric"))

# Line search algorithm from Hager & Zhang ???
# Adapted from LineSearch.jl
HagerZhang <- function (dphifn, phi_0, dphi_0, control = HagerZhangControl())
{
  nextfloat <- function(x)
    x + .Machine$double.eps

  for (var in names(control)) assign(var, control[[var]])

  if (!(is.finite(phi_0) && is.finite(dphi_0)))
    stop("Value and slope at step length = 0 must be finite.")
  if (dphi_0 >= .Machine$double.eps * abs(phi_0))
    stop("Search direction is not a direction of descent.")
  else if (dphi_0 >= 0)
    return(0) #return(list(0, phi_0))

  # Prevent values of x_new = x+alpha*s that are likely to make
  # phi(x_new) infinite
  iterfinitemax = ceiling(-log(.Machine$double.eps, 2))
  st <- new("HagerZhangStorage", alphas = c(0), values = c(phi_0), slopes = c(dphi_0))
  if (verbose)
    cat("New linesearch\n")

  phi_lim = phi_0 + epsilon * abs(phi_0)
  stopifnot(c >= 0)
  c <= .Machine$double.eps && return(0) #return(list(0, phi_0))
  stopifnot(is.finite(c) && c <= alphamax)
  fit = dphifn(c)
  phi_c = fit$objective
  dphi_c = fit$gradient
  iterfinite = 1
  while(!(is.finite(phi_c) && is.finite(dphi_c)) && iterfinite < iterfinitemax)
  {
    #mayterminate = FALSE #bc c is set to 1 for now
    iterfinite = iterfinite + 1
    c = c * psi3
    fit = dphifn(c)
    phi_c = fit$objective
    dphi_c = fit$gradient
  }
  if (!(is.finite(phi_c) && is.finite(dphi_c)))
  {
    warning("Failed to achieve finite new evaluation point, using alpha=0")
    return(0) #return(list(0, phi0))
  }
  st$alphas <- c(st$alphas, c)
  st$values <- c(st$values, phi_c)
  st$slopes <- c(st$slopes, dphi_c)

  #NSP: use initial c = 1 for Newton method. Uncommenting below slows things down
  # If c was generated by quadratic interpolation, check whether it
  # satisfies the Wolfe conditions
  #if (mayterminate && satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma))
  #if (satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma))
  #{
  #  if(verbose)
  #    cat("Wolfe condition satisfied on point alpha = ", c, "\n")
  #  return(list(c, phi_c)) # phi_c
  #}

  # Initial bracketing step (HZ, stages B0-B3)
  isbracketed = FALSE
  ia = 1
  ib = 2
  stopifnot(length(st$alphas) == 2)
  iter = 1
  cold = -1
  while(!isbracketed && iter < linesearchmax)
  {
    if(verbose)
      cat("bracketing: ia = ", ia, ", ib = ", ib, ", c = ", c, ", phi_c = ", phi_c, ", dphi_c = ", dphi_c, "\n")
    if(dphi_c >= 0)
    {
      # We've reached the upward slope, so we have b; examine
      # previous values to find a
      ib = length(st$alphas)
      for (i in (ib - 1):1)
      {
        if (st$values[i] <= phi_lim)
        {
          ia = i
          break
        }
      }
      isbracketed = TRUE
    }
    else if (st$values[length(st$values)] > phi_lim)
    {
      # The value is higher, but the slope is downward, so we must
      # have crested over the peak. Use bisection.
      ib = length(st$alphas)
      ia = 1
      if (c !=  st$alphas[ib] || st$slopes[ib] >= 0)
        stop("c = ", c)
      #bis = bisect(dphifn, alphas, values, slopes, ia, ib, phi_lim, verbose)
      bis = bisect(dphifn, st, ia, ib, phi_lim, verbose)
      ia  = bis[[1]]
      ib  = bis[[2]]
      isbracketed = TRUE
    }
    else
    {
      # We'll still going downhill, expand the interval and try again.
      # Reaching this branch means that dphi_c < 0 and phi_c <= phi_0 + e_k 
      # So cold = c has a lower objective than phi_0 up to epsilon. 
      # This makes it a viable step to return if bracketing fails.

      # Bracketing can fail if no cold < c <= alphamax can be found with finite phi_c and dphi_c. 
      # Going back to the loop with c = cold will only result in infinite cycling.
      # So returning (cold, phi_cold) and exiting the line search is the best move.
      cold = c
      phi_cold = phi_c
      if (nextfloat(cold) >= alphamax)
      {
        return(cold) #list(c(cold, phi_cold))
      }
      c = c * rho
      if (c > alphamax)
      {
        c = alphamax
        if (verbose)
          cat("bracket: exceeding alphamax, using c = alphamax = ", alphamax, ", cold = ", cold, "\n")
      }
      fit = dphifn(c)
      phi_c = fit$objective
      dphi_c = fit$gradient
      iterfinite = 1
      while (!(is.finite(phi_c) && is.finite(dphi_c)) && c > nextfloat(cold) && iterfinite < iterfinitemax)
      {
        alphamax = c # shrinks alphamax, assumes that steps >= c can never have finite phi_c and dphi_c
        iterfinite = iterfinite + 1
        if (verbose)
          cat("bracket: non-finite value, bisection\n")
        c = (cold + c) / 2
        fit = dphifn(c)
        phi_c = fit$objective
        dphi_c = fit$gradient
      }
      if(!(is.finite(phi_c) && is.finite(dphi_c)))
      {
        if (verbose)
        {
          cat("Warning: failed to expand interval to bracket with finite values. If this happens frequently, check your function and gradient.\n")
          cat("c = ", c, ", alphamax = ", alphamax, ", phi_c = ", phi_c, ", dphi_c = ", dphi_c, "\n")
        }
        return(cold) #return(list(cold, phi_cold))
      }
      st$alphas = c(st$alphas, c)
      st$values = c(st$values, phi_c)
      st$slopes = c(st$slopes, dphi_c)
    }
    iter = iter + 1
  }
  while (iter < linesearchmax)
  {
    a = st$alphas[ia]
    b = st$alphas[ib]
    stopifnot(b > a)
    if (verbose)
      cat("linesearch: ia = ", ia, ", ib = ", ib, ", a = ", a, ", b = ", b, ", phi(a) = ", st$values[ia], ", phi(b) = ", st$values[ib], "\n")
    if (b - a <= .Machine$double.eps)
    {
      return(a) #return(c(a, st$values[ia])) 
    }
    #sec = secant2(dphifn, alphas, values, slopes, ia, ib, phi_lim, delta, sigma, verbose)
    sec = secant2(dphifn, st, ia, ib, phi_lim, delta, sigma, verbose)
    iswolfe = sec[[1]]
    iA = sec[[2]]
    iB = sec[[3]]
    if(iswolfe)
    {
      return (st$alphas[iA]) #return(list(st$alphas[iA], st$values[iA]))
    }
    A = st$alphas[iA]
    B = st$alphas[iB]
    stopifnot(B > A)
    if (B - A < gamma * (b - a))
    {
      if (verbose)
        cat("Linesearch: secant succeeded\n")
      if (nextfloat(st$values[ia]) >= st$values[ib] && nextfloat(st$values[iA]) >= st$values[iB])
      {
        # It's so flat, secant didn't do anything useful, time to quit
        if (verbose)
          cat("Linesearch: secant suggests it's flat\n")
        return(A) #return(list(A, st$values[iA]))
      }
      ia = iA
      ib = iB
    }
    else
    {
      # Secant is converging too slowly, use bisection
      if(verbose)
        cat("Linesearch: secant failed, using bisection\n")
      c = (A + B) / 2

      fit = dphifn(c)
      phi_c = fit$objective
      dphi_c = fit$gradient
#      while (!(is.finite(phi_c) && is.finite(dphi_c))) #NSP: sometimes midpoint is not OK even if either end of interval is
#      {
#        c = (A + c) / 2
#        fit = dphifn(c)
#        phi_c = fit$objective
#        dphi_c = fit$gradient
#      } #end NSP
      stopifnot(is.finite(phi_c) && is.finite(dphi_c))
      st$alphas = c(st$alphas, c)
      st$values = c(st$values, phi_c)
      st$slopes = c(st$slopes, dphi_c)

      upd = update(dphifn, st, iA, iB, length(st$alphas), phi_lim, verbose) 
      ia = upd[[1]]
      ib = upd[[2]]
    }
    iter = iter + 1
  }
  stop("Linesearch failed to converge, reached maximum iterations ", st$alphas[ia])
}

# Check Wolfe & approximate Wolfe
satisfies_wolfe <- function(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma)
{
    wolfe1 = delta * dphi_0 >= (phi_c - phi_0) / c && dphi_c >= sigma * dphi_0
    wolfe2 = (2 * delta - 1) * dphi_0 >= dphi_c && dphi_c >= sigma * dphi_0 && phi_c <= phi_lim
    return (wolfe1 || wolfe2)
}

# HZ, stages S1-S4
secant0 <- function(a, b, dphi_a, dphi_b)
{
  return ((a * dphi_b - b * dphi_a) / (dphi_b - dphi_a))
}

secant <- function(alphas, values, slopes, ia, ib) #erg, remove this and insert indexing
{
  return (secant0(alphas[ia], alphas[ib], slopes[ia], slopes[ib]))
}

# phi
#secant2 <- function(dphifn, alphas, values, slopes, ia, ib, phi_lim, delta, sigma, verbose)
secant2 <- function(dphifn, st, ia, ib, phi_lim, delta, sigma, verbose)
{
    phi_0 = st$values[1]
    dphi_0 = st$slopes[1]
    a = st$alphas[ia]
    b = st$alphas[ib]
    dphi_a = st$slopes[ia]
    dphi_b = st$slopes[ib]
    if(!(dphi_a < 0 && dphi_b >= 0))
        stop("Search direction is not a direction of descent; this error may indicate that user-provided derivatives are inaccurate.\n(dphi_a = ", dphi_a, "dphi_b = ", dphi_b)
    c = secant0(a, b, dphi_a, dphi_b)
    if (verbose)
      cat("secant2: a = ", a, ", b = ", b, ", c = ", c, "\n")
    stopifnot(is.finite(c))
    fit = dphifn(c)
    phi_c = fit$objective
    dphi_c = fit$gradient
    stopifnot(is.finite(phi_c) && is.finite(dphi_c))

    st$alphas = c(st$alphas, c)
    st$values = c(st$values, phi_c)
    st$slopes = c(st$slopes, dphi_c)

    ic = length(st$alphas)
    if (satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma))
    {
      if (verbose)
        cat("secant2: first c satisfied Wolfe conditions\n")
      return(list(TRUE, ic, ic))
    }

    upd = update(dphifn, st, ia, ib, ic, phi_lim, verbose) 
    iA = upd[[1]]
    iB = upd[[2]]
    if (verbose)
      cat("secant2: iA = ", iA, ", iB = ", iB, ", ic = ", ic, "\n")
    a = st$alphas[iA]
    b = st$alphas[iB]
    doupdate = FALSE
    if (iB == ic)
    {
      # we updated b, make sure we also update a
      c = secant(st$alphas, st$values, st$slopes, ib, iB)
    }
    else if (iA == ic)
    {
      # we updated a, do it for b too
      c = secant(st$alphas, st$values, st$slopes, ia, iA)
    }
    if ((iA == ic || iB == ic) && a <= c && c <= b)
    {
        if (verbose)
            cat("secant2: second c = ", c, "\n")
        fit = dphifn(c)
        phi_c = fit$objective
        dphi_c = fit$gradient
        stopifnot(is.finite(phi_c) && is.finite(dphi_c))

        st$alphas = c(st$alphas, c)
        st$values = c(st$values, phi_c)
        st$slopes = c(st$slopes, dphi_c)

        ic = length(st$alphas)
        # Check arguments here
        if (satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma))
        {
            if (verbose)
                cat("secant2: second c satisfied Wolfe conditions\n")
            return (list(TRUE, ic, ic))
        }
        upd = update(dphifn, st, iA, iB, ic, phi_lim, verbose) 
        iA = upd[[1]]
        iB = upd[[2]]
    }
    if (verbose)
        cat("secant2 output: a = ", st$alphas[iA], ", b = ", st$alphas[iB], "\n")
    return(list(FALSE, iA, iB))
}

# HZ, stages U0-U3
# Given a third point, pick the best two that retain the bracket
# around the minimum (as defined by HZ, eq. 29)
# b will be the upper bound, and a the lower bound
#update <- function(dphifn, alphas, values, slopes, ia, ib, ic, phi_lim, verbose)
update <- function(dphifn, st, ia, ib, ic, phi_lim, verbose)
{
    a = st$alphas[ia]
    b = st$alphas[ib]
    # Debugging (HZ, eq. 4.4):
    stopifnot(st$slopes[ia] < 0)
    stopifnot(st$values[ia] <= phi_lim)
    stopifnot(st$slopes[ib] >= 0)
    stopifnot(b > a)
    c = st$alphas[ic]
    phi_c = st$values[ic]
    dphi_c = st$slopes[ic]
    if (verbose)
        cat("update: ia = ", ia, ", a = ", a, ", ib = ", ib, ", b = ", b, ", c = ", c, ", phi_c = ", phi_c, ", dphi_c = ", dphi_c)
    if (c < a || c > b)
        return (list(ia, ib)) #, 0, 0  # it's out of the bracketing interval
    if (dphi_c >= 0)
        return (list(ia, ic)) #, 0, 0  # replace b with a closer point
    # We know dphi_c < 0. However, phi may not be monotonic between a
    # and c, so check that the value is also smaller than phi_0.  (It's
    # more dangerous to replace a than b, since we're leaving the
    # secure environment of alpha=0; that's why we didn't check this
    # above.)
    if (phi_c <= phi_lim)
        return (list(ic, ib))#, 0, 0  # replace a
    # phi_c is bigger than phi_0, which implies that the minimum
    # lies between a and c. Find it via bisection.
    #return (bisect(dphifn, alphas, values, slopes, ia, ic, phi_lim, verbose))
    return (bisect(dphifn, st, ia, ic, phi_lim, verbose))
}

# HZ, stage U3 (with theta=0.5)
#bisect <- function(dphifn, alphas, values, slopes, ia, ib, phi_lim, verbose)
bisect <- function(dphifn, st, ia, ib, phi_lim, verbose)
{
    gphi = NaN
    a = st$alphas[ia]
    b = st$alphas[ib]
    # Debugging (HZ, conditions shown following U3)
    stopifnot(st$slopes[ia] < 0)
    stopifnot(st$values[ia] <= phi_lim)
    stopifnot(st$slopes[ib] < 0)       # otherwise we wouldn't be here
    stopifnot(st$values[ib] > phi_lim)
    stopifnot(b > a)
    while (b - a > .Machine$double.eps)
    {
        if (verbose)
            cat("bisect: a = ", a, ", b = ", b, ", b - a = ", b - a, "\n")
        d = (a + b) / 2.
        fit = dphifn(d)
        phi_d = fit$objective
        gphi = fit$gradient
#        if (!(is.finite(phi_d) && is.finite(gphi))) #NSP: sometimes we have that function is well-behaved at endpoints but not intermediate
#        {
#          b = d
#          next
#        } #end NSP
        stopifnot(is.finite(phi_d) && is.finite(gphi))

        st$alphas = c(st$alphas, d)
        st$values = c(st$values, phi_d)
        st$slopes = c(st$slopes, gphi)

        id = length(st$alphas)
        if (gphi >= 0)
        {
            return (list(ia, id)) # replace b, return
        }
        if (phi_d <= phi_lim)
        {
            a = d # replace a, but keep bisecting until dphi_b > 0
            ia = id
        }
        else
        {
            b = d
            ib = id
        }
    }
    return (list(ia, ib))
}

#' Control settings for Newton-Raphson optimizer
#'
#' TODO after finalizing line search
#'
#' @export
NewtonRaphsonControl <- function(maxit = 100, ctol = sqrt(.Machine$double.eps), ftol = sqrt(.Machine$double.eps), etol = 10*.Machine$double.eps, verbose = FALSE, eps = 1e-8, ls.control = HagerZhangControl())
  list(maxit = maxit, ctol = ctol, etol = etol, ftol = ftol, verbose = verbose, eps = eps, ls.control = ls.control)

# deprecated, should delete at some point
NewtonRaphson <- function(par, fn, control = NewtonRaphsonControl())
{
  NewtonRaphsonNaN <- function()
  {
    list(objective = NaN,
         gradient  = matrix(NaN, length(par), 1),
         hessian   = matrix(NaN, length(par), length(par)))
  }

  for (var in names(control)) assign(var, control[[var]])
  etol <- etol * length(par)

  if (verbose)
    cat("Newton-Raphson with Hager-Zhang line search\n")

  convergence <- 0
  oldfit      <- list(objective = Inf)
  oldpar      <- par

  for (i in 1:maxit)
  {
    fit <- fn(par, gradient = TRUE, hessian = TRUE)

    if (verbose)
      cat(paste0("[", i, "]"), 
          "f(x) =", formatC(-fit$objective, digits=8, width=8, format="f"), 
          "  delta f(x) =", if(i > 1) formatC(oldfit$objective - fit$objective, digits=8, width=8, format="f") else NaN, 
          "  max |f'(x)| =", formatC(max(abs(fit$gradient)), digits=8, width=8, format="f"),
          "  det -f''(x) =", formatC(det(fit$hessian), digits=8, width=8, format="f"),
          "\n")

    if (max(abs(fit$gradient)) < ctol || (i > 1 && oldfit$objective - fit$objective < ftol))
      break

    ehess        <- eigen(fit$hessian)
    ehess$values <- abs(ehess$values)
    ehess$values <- ifelse(ehess$values < max(abs(fit$hessian)) * etol, 1, ehess$values)
    ihess        <- ehess$vectors %*% solve(diag(ehess$values, nrow=length(par))) %*% t(ehess$vectors)
    des          <- -ihess %*% fit$gradient
    phi0         <- fit$objective
    dphi0        <- c(t(des) %*% fit$gradient)

    dphi_fn <- function(alpha) 
    {
      tryCatch({
        phi <- fn(par + alpha*des, gradient = TRUE, hessian = FALSE)
        list(objective = phi$objective, gradient = c(t(des) %*% phi$gradient))
      }, error = function(e) {
        NewtonRaphsonNaN()
      })
    }

    alpha <- HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    par   <- par + des * alpha

    oldfit       <- fit
  }
  if (verbose)
    cat ("Finished with max(abs(grad)) =", max(abs(fit$gradient)), "and delta f(x) =", oldfit$objective - fit$objective, "\n")

  if (i == maxit)
  {
    warning("Maxit reached for Newton steps")
    convergence = 1
  } 

  list(par = par,
       gradient = fit$gradient,
       hessian = fit$hessian,
       value = fit$objective,
       fit = fit,
       convergence = convergence)
}

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

#' Fst from allele counts
#'
#' Calculates Fst using the estimator of Bhatia et al (2015)
#' from counts of the derived allele across biallelic markers
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @return A matrix containing pairwise Fst
#'
#' @references
#' Bhatia et al 2015 TODO
#'
#' @export

fst_from_biallelic <- function(Y, N)
{
  # TODO: missing data

  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")

  Fr <- Y/N
  f2 <- apply(apply(Fr, 2, function(x) outer(x,x,"-")^2), 1, mean, na.rm=TRUE)
  cornum <- apply(Fr*(1-Fr)/(N-1), 1, mean, na.rm=TRUE)
  corden <- apply(Fr*(1-Fr)*N/(N-1), 1, mean, na.rm=TRUE)
  fst <- (f2 - c(outer(cornum, cornum, "+")))/(f2 - c(outer(cornum, cornum, "+")) + c(outer(corden, corden, "+")))
  fst <- matrix(fst, nrow(Y), nrow(Y))
  diag(fst) <- 0
  fst
}

#' Covariance from allele counts
#'
#' Calculates covariance from counts of the derived allele across biallelic markers,
#' using normalized allele frequencies
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @details
#' TODO formula
#'
#' @return A matrix containing pairwise covariance
#'
#' @export

cov_from_biallelic <- function(Y, N)
{
  #TODO missing data

  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")

  Fr <- colSums(Y) / colSums(N)
  Y  <- (Y - N*Fr) / sqrt(N * Fr * (1-Fr))

  Y %*% t(Y) / ncol(Y)
}

#' Distance matrix from covariance matrix
#'
#' Returns the squared-distance matrix associated with a given covariance matrix
#'
#' @param Cov A covariance matrix (does not need to be full-rank)
#'
#' @details
#' TODO formula
#'
#' @return A distance matrix of the same dimension as the input
#'
#' @export

dist_from_cov <- function(Cov)
{
  ones <- matrix(1, nrow(Cov), 1)
  diag(Cov) %*% t(ones) + ones %*% t(diag(Cov)) - 2*Cov
}

#' Distance from allele counts
#'
#' Calculates genetic distance from counts of the derived allele across biallelic markers,
#' using normalized allele frequencies
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @details
#' TODO formula
#'
#' @return A matrix containing pairwise distance
#'
#' @export

dist_from_biallelic <- function(Y, N)
{
  dist_from_cov(cov_from_biallelic(Y, N))
}

# TODO check and document
wishart_simulate_distance <- function(seed, S, nu)
{
  set.seed(seed)
  dist_from_cov(solve(rWishart(1, nu, S)[,,1]))
}

wishart_simulate_experiment <- function(seed, N, P, K, nu, neval, timingOnly=FALSE)
{
  library(nloptr)
  library(RandomFields)

  covariates <- list()
  for(k in 1:K)
  {
    covariates[[paste0("var", k)]] <- raster::raster(matrix(RandomFields::RFsimulate(RandomFields::RMexp(scale=10), expand.grid(x=1:N, y=1:N), seed=seed)@data[[1]],N,N))
    covariates[[paste0("var", k)]] <- covariates[[paste0("var", k)]]/max(covariates[[paste0("var", k)]][])
  }
  covariates <- raster::stack(covariates)

  set.seed(seed)
  demes <- sample(1:(N^2), P)
  coords <- raster::xyFromCell(covariates[[1]], demes, spatial = TRUE)

  surf <- radish_conductance_surface(covariates, coords, directions = 4, saveStack = FALSE)

  beta <- matrix(rnorm(K), 1, K)
  E <- radish_distance(loglinear_conductance, surf, beta, covariance = TRUE)$covariance[,,1]
  S <- wishart_simulate_distance(seed=seed, nu=nu, S=E)
  S <- S/(max(S)*10)

  # functions
  nograd <- function(par)
  {
    radish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = FALSE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  wgrad <- function(par)
  {
    radish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  whess <- function(par)
  {
    radish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = TRUE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }

  timings_whess <- c()
  timings_wgrad <- c()
  timings_nograd <- c()
  for(i in 1:neval)
  {
    timings_whess <- rbind(timings_whess, system.time(ll <- whess(rep(0,K))))
    timings_wgrad <- rbind(timings_wgrad, system.time(ll <- wgrad(rep(0,K))))
    timings_nograd <- rbind(timings_nograd, system.time(ll <- nograd(rep(0,K))))
  }

  #optimization
  if(!timingOnly){
  opt_newton <- radish_optimize(loglinear_conductance, leastsquares, surf, S, leverage = FALSE)

  .fcall <<- 0 
  bqfn <- function(par)
  {
    .fcall <<- .fcall + 1 #debug
    nograd(par)$objective
  }
  opt_bobyqa <- bobyqa(rep(0, K), bqfn)
  opt_bobyqa$fcall <- .fcall
} else {
  opt_bobyqa = list(fcall=NA, value=NA)
  opt_newton = list(fcall=NA, loglik=NA, fit=list(boundary=NA))
  }


  list(timings=list(nograd=timings_nograd, wgrad=timings_wgrad, whess=timings_whess), opt_newton=opt_newton, 
       #opt_lbfgs=opt_lbfgs, 
       opt_bobyqa=opt_bobyqa, seed=seed, beta=beta, K=K, N=N, P=P)
}

run_benchmarks <- function(K = c(1,2,4,8,16,32), N = c(100), P = c(30), reps=5, timingOnly=timingOnly)
{
  set.seed(1)
  seeds <- sample.int(100000, length(K)*length(N)*length(P)*reps)
  z <- 0
  out <- list()
  for(p in P)
  {
    out[[as.character(p)]] <- list()
    for(n in N)
    {
      out[[as.character(p)]][[as.character(n)]] <- list()
      for(k in K)
      {
        out[[as.character(p)]][[as.character(n)]][[as.character(k)]] <- list()
        for(rep in 1:reps)
        {
          z <- z + 1
          try({
            out[[as.character(p)]][[as.character(n)]][[as.character(k)]][[as.character(rep)]] <- 
              wishart_simulate_experiment(seeds[z], n, p, k, 100, 1, timingOnly=timingOnly)
          })
        }
      }
    }
  }
  out
}
#aahh <- run_benchmarks(c(1,2,4,8,16,32), c(200), c(30), 10, timingOnly=TRUE)

extract_benchmark_timing <- function(benchmarks)
{
  out <- c()
  for(p in names(benchmarks))
    for(n in names(benchmarks[[p]]))
      for(k in names(benchmarks[[p]][[n]]))
        for(r in names(benchmarks[[p]][[n]][[k]]))
        {
          tmp <- benchmarks[[p]][[n]][[k]][[r]]
          oo <- data.frame(P=as.numeric(p), N=as.numeric(n), K=as.numeric(k), rep=as.numeric(r),
                           timing_nograd = mean(tmp$timings$nograd[,1]),
                           timing_wgrad = mean(tmp$timings$wgrad[,1]),
                           timing_whess = mean(tmp$timings$whess[,1]),
                           newton_eval  = tmp$opt_newton$fcall,
                           newton_ll    = tmp$opt_newton$loglik,
                           newton_boundary = tmp$opt_newton$fit$boundary,
                           #lbfgs_eval  = tmp$opt_lbfgs$fcall,
                           bobyqa_eval  = tmp$opt_bobyqa$fcall,
                           bobyqa_ll    = tmp$opt_bobyqa$value
                           )
          out <- rbind(out, oo)
        }
  out
}

#' Create a parameterized conductance surface
#'
#' Given a set of spatial covariates and a set of spatial coordinates, create a
#' graph representing a parameterized conductance surface.
#'
#' @param covariates A \code{RasterStack} containing spatial covariates
#' @param coords A \code{SpatialPoints} object containing coordinates for a set of focal cells, with the same projection as the rasters in \code{covariates}
#' @param directions If \code{4}, consider only horizontal/vertical neighbours as adjacent; if \code{8}, also consider diagonal neighbours as adjacent
#' @param saveStack If \code{TRUE}, the \code{RasterStack} is returned with missing data masked uniformly across rasters
#'
#' @details NAs are shared across rasters in \code{covariates}, and a warning is thrown if a given cell has mixed NA and non-NA values across the stack. Comparing models with different patterns of missing spatial data (e.g. fit to different stacks of rasters) can give superficially inconsistant results, as these essentially involve different sets of vertices. Thus model comparison should use models fitted to the same \code{radish_graph} object.
#'
#' Disconnected components are identified and removed, so that only the largest connected component in the graph is retained. The function aborts if there are focal cells that belong to a disconnected component.
#'
#' Rasters of categorical covariates must have an associated RAT to be 'recognized' as categorical, see \code{\link[raster]{ratify}} and \code{examples} below. The names of levels are taken from the \code{VALUE} column of the RAT, if it exists (otherwise, the integer codes in the \code{ID} column are used).
#'
#' @seealso \code{\link{radish}}, \code{\link[raster]{stack}}
#'
#' @return An object of class \code{radish_graph}
#'
#' @references
#'
#' Pope NS. In prep. Fast gradient-based optimization of resistance surfaces.
#'
#' @examples
#'
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=raster::scale(melip.altitude), 
#'                                  forestcover=raster::scale(melip.forestcover)))
#'
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' # categorical covariates:
#' # rasters of categorical covariates must have an associated RAT, see 'details'
#' forestcover_class <- cut(raster::values(melip.forestcover), breaks = c(0, 1/3, 1/6, 1)) 
#' melip.forestcover_cat <- 
#'   raster::ratify(raster::setValues(melip.forestcover, as.numeric(forestcover_class)))
#'
#' RAT <- levels(melip.forestcover_cat)[[1]]
#' RAT$VALUE <- levels(forestcover_class) #explicitly define level names
#' levels(melip.forestcover_cat) <- RAT
#' 
#' covariates_cat <- raster::stack(list(forestcover = melip.forestcover_cat,
#'                                      altitude = melip.altitude)) 
#'
#' @export

conductance_surface <- function(covariates, coords, directions=4, saveStack=TRUE)
{
  stopifnot(class(covariates) == "RasterStack")
  stopifnot(class(coords) == "SpatialPoints")
  stopifnot(directions %in% c(4, 8))

  # throw a warning if missing cells are not identical across layers
  spdat   <- as.matrix(raster::getValues(covariates))
  missing <- is.na(rowSums(spdat))
  if (!all(apply(spdat[missing,,drop=FALSE], 1, function(x) all(is.na(x)))))
    warning("Missing cells are not identical across rasters; be careful regarding model selection (see ?conductance_surface)")

  # share missing cells across layers and ensure graph is fully connected
  if (any(missing))
  {
    spdat[missing,] <- NA
    cr <- covariates[[1]]
    raster::values(cr) <- ifelse(missing, NA, 1) #b/c of how clump handles zeros
    cr <- raster::clump(cr, directions = directions)
    connected_component <- names(which.max(table(raster::getValues(cr)))) #using the largest connected subgraph
    disconnected <- raster::getValues(cr) != as.integer(connected_component)
    disconnected[is.na(disconnected)] <- FALSE
    if (any(disconnected))
      warning(paste("Pruned", sum(disconnected), "disconnected cells across rasters"))
    missing <- disconnected | missing
    spdat[missing,] <- NA
    for (i in 1:dim(covariates)[3])
      raster::values(covariates[[i]]) <- spdat[,i]
  }

  # get adjacency list
  adj <- raster::adjacent(covariates[[1]], which(!missing), target=which(!missing), directions=directions)

  # find cells of demes
  cells <- unmapped_cells <- raster::cellFromXY(covariates[[1]], coords)

  # remove NAs and remap indices of adjacency list to be contiguous
  map     <- cbind(1:length(which(!missing)), (1:raster::ncell(covariates[[1]]))[which(!missing)])
  spdat   <- spdat[!missing,,drop=FALSE]
  adj[,1] <- map[match(adj[,1], map[,2]),1]
  adj[,2] <- map[match(adj[,2], map[,2]),1]
  cells   <- map[match(cells, map[,2]),1]

  # check that cells lie on connected portion of raster
  if(any(is.na(raster::getValues(covariates[[1]])[unmapped_cells])))
    stop("At least one deme is located on a missing cell")

  # figure out which raster layers are factors
  is_factor <- sapply(names(covariates), function(i) raster::is.factor(covariates[[i]]))
  spdat <- as.data.frame(spdat)
  if (any(is_factor))
  {
    factors <- names(covariates)[is_factor]
    warning("Treating covariates \"", paste(factors, collapse = "\" \""), "\" as factors")
    for (i in factors)
    {
      levels    <- raster::levels(covariates[[i]])[[1]]
      ids       <- levels$ID
      levels    <- if(is.null(levels$VALUE)) ids else levels$VALUE
      spdat[,i] <- factor(levels[match(spdat[,i], ids)])
    }
  }

  # form and factorize Laplacian
  N    <- nrow(spdat)
  Q    <- Matrix::sparseMatrix(i = adj[,1], j = adj[,2], dims = c(N, N),
                               x = -rep(1, nrow(adj)),
                               use.last.ij = TRUE)
  Q    <- Matrix::forceSymmetric(Q)
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  adj  <- rbind(Q@i, rep(1:Q@Dim[2] - 1, diff(Q@p))) #upper-triangular, 0-based
  LQn  <- Matrix::Cholesky(Qn, LDL = TRUE) 

  out <- list("demes"        = cells,
              "x"            = spdat,
              "adj"          = adj,
              "covariates"   = colnames(spdat), 
              "laplacian"    = Q,
              "choleski"     = LQn,
              "stack"        = if(saveStack) covariates else NULL)
  class(out) <- "radish_graph"
  out
}

conductance <- function(x, ...)
{
  UseMethod("conductance")
}

conductance.radish_graph <- function(object, fit, quantile = 0.95)
{
  stopifnot(class(fit) == "radish")
  stopifnot(!fit$fit$boundary && !is.null(fit$mle$theta))

  conductance <- fit$submodels$f(fit$mle$theta)
  ci <- conductance$confint(theta = fit$mle$theta, vcov = -solve(fit$mle$hessian), 
                            quantile = quantile, scale = "conductance")
  colnames(ci) <- paste0(c("lower", "upper"), round(100*quantile, 1))

  if (!is.null(object$stack))
  {
    template <- object$stack[[1]]
    missing  <- is.na(raster::values(template))
    raster::values(template)[!missing] <- 1.
    lower <- upper <- est <- template
    raster::values(est)[!missing] <- conductance$conductance
    raster::values(lower)[!missing] <- ci[,1]
    raster::values(upper)[!missing] <- ci[,2]
    out <- list()
    out[["est"]] <- est
    out[[colnames(ci)[1]]] <- lower
    out[[colnames(ci)[2]]] <- upper
    raster::stack(out)
  }
  else
  {
    warning("No rasters associated with graph, returning conductance as vector")
    out <- cbind("est" = conductance$conductance, ci)
    out
  }
}

#for validation and debugging, generate a simple 1D lattice
Lattice1D <- function(spdat, coords, fn)
{
  stopifnot(all(!is.na(spdat)))
  stopifnot(all(coordinates(coords) >= 0 & coordinates(coords) <= 1))

  rl <- list()
  for(i in 1:ncol(spdat))
    rl[[paste0("var", i)]] <- raster::raster(as.matrix(spdat[,i]))
  rl <- raster::stack(rl)

  conductance_surface(rl, coords)
}

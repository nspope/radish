#' Create a parameterized conductance surface
#'
#' Given a set of spatial covariates (as a 'RasterStack') and a set of
#' spatial coordinates (class 'SpatialPoints'), create an abstract parameterized
#' conductance surface.
#'
#' @param formula A 'formula' object for the spatial data in 'covariates' (see 'Details')
#' @param covariates A 'RasterStack' containing spatial covariates
#' @param coords A 'SpatialPoints' object containing coordinates for a set of focal cells
#' @param directions If '4', consider only horizontal/vertical neighbours as adjacent; if '8', also consider diagonal neighbours as adjacent
#' @param saveStack If TRUE, the modified RasterStack is returned
#'
#' @details NAs are shared across rasters in 'covariates', and a warning is thrown if a given cell has mixed NA and non-NA values across the stack. Comparing models with different patterns of missing spatial data (e.g. fit to different stacks of rasters) can give superficially inconsistant results, as these essentially involve different sets of vertices. The recommended way to create sets of nested models (with an invariant set of vertices) is via '?downdate'.
#'
#' Disconnected components are identified and removed, so that only the largest connected component in the graph is retained. The function aborts if there are focal cells that belong to a disconnected component.
#'
#' The design matrix for the spatial data that is ultimately used in 'radish' (by an object of class 'conductance_model') is created from 'formula' and 'covariates' via 'model.matrix'.
#'
#' @return An object of class 'radish_graph'
#'
#' @references
#'
#' Pope NS. In prep. Fast gradient-based optimization of resistance surfaces.
#'
#' @examples
#' library(raster)
#' 
#' data(melip)
#' 
#' covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
#' surface <- radish_conductance_surface(~altitude * forestcover, covariates, melip.coords, directions = 8)
#'
#' @export

conductance_surface <- function(formula, covariates, coords, directions=4, saveStack=TRUE)
{
  stopifnot(class(formula) == "formula")
  stopifnot(class(covariates) == "RasterStack")
  stopifnot(class(coords) == "SpatialPoints")
  stopifnot(directions %in% c(4, 8))

  # check if formula is consistant with data, remove response, add intercept
  formula_covariates <- attr(delete.response(terms(formula)), "factors")
  stopifnot(rownames(formula_covariates) %in% names(covariates))
  formula <- reformulate(colnames(formula_covariates))

  # if any layers are not in formula, remove them
  missing_covariates <- !(names(covariates) %in% rownames(formula_covariates))
  if (any(missing_covariates))
  {
    unused_covariates <- names(covariates)[missing_covariates]
    warning("Removed unused spatial covariates: ", 
            paste(unused_covariates, collapse = " "))
    covariates <- raster::dropLayer(covariates, unused_covariates)
  }

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

  # get model matrix and check for rank deficiency
  # NOTE: in future I could make this sparse via Matrix::sparse.model.matrix?
  spdat <- model.matrix(formula, data = spdat)
  stopifnot(qr(spdat)$rank == ncol(spdat))
  spdat <- spdat[,colnames(spdat) != "(Intercept)", drop=FALSE]

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
              "formula"      = formula,
              "stack"        = if(saveStack) covariates else NULL)
  class(out) <- "radish_graph"
  out
}

downdate <- function(x, ...)
{
  UseMethod("downdate")
}

downdate.radish_graph <- function(object, formula) 
{
  # update formula
  object$formula <- update(object$formula, formula)

  # update model matrix
  object$x <- model.matrix(object$formula, data = data.frame(object$x))
  object$x <- object$x[,colnames(object$x) != "(Intercept)",drop=FALSE]

  object$covariates <- colnames(object$x)

  object
} #NOTE: this can only 'downdate' an existing 'radish_graph'

conductance <- function(x, ...)
{
  UseMethod("conductance")
}

conductance.radish_graph <- function(object, fit, ci = 0.95)
{
  stopifnot(class(fit) == "radish")
  stopifnot(!fit$boundary)

  if (!is.null(object$stack))
  {
    template <- object$stack[[1]]
    missing  <- is.na(values(template))
    values(template)[!missing] <- 1.
    # TODO: confidence intervals
    #lower <- upper <- est <- template
    est <- template
    values(est)[!missing] <- fit$f(object$x, fit$theta)$conductance
    return(est)
  }
  else
  {
    warning("No raster associated with graph, returning conductance as vector")
    return(fit$f(object$x, fit$theta)$conductance)
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

  conductance_surface(reformulate(names(rl)), rl, coords)
}

#' Create a parameterized conductance surface
#'
#' Given a set of spatial covariates (as a 'RasterStack') and a set of
#' spatial coordinates (class 'SpatialPoints'), create an abstract parameterized
#' conductance surface.
#'
#' @param covariates A 'RasterStack' containing spatial covariates
#' @param coords A 'SpatialPoints' object containing coordinates for a set of target cells
#' @param directions If '4', consider only horizontal/vertical neighbours as adjacent; if '8', also consider diagonal neighbours as adjacent
#' @param saveStack If TRUE, the modified RasterStack is returned
#'
#' @details NAs are shared across components of the RasterStack, and a warning is thrown if a given cell has mixed NA and non-NA values across the stack. Disconnect components are identified and removed (only the largest connected component is retained).
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
#' surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
#'
#' @export

radish_conductance_surface <- function(covariates, coords, directions=4, saveStack=TRUE)
{

  stopifnot(class(covariates) == "RasterStack"  )
  stopifnot(    class(coords) == "SpatialPoints")
  stopifnot(directions %in% c(4, 8))

  # throw a warning if missing cells are not identical across layers
  spdat   <- as.matrix(raster::getValues(covariates))
  missing <- is.na(rowSums(spdat))
  if (!all(apply(spdat[missing,,drop=FALSE], 1, function(x) all(is.na(x)))))
    warning("Missing cells are not identical across rasters: model selection may be untrustworthy")

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
  LQn  <- Matrix::Cholesky(Qn, LDL = TRUE) #TODO: supernodal?

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

#for validation and debugging, generate a simple 1D lattice
Lattice1D <- function(spdat, coords, fn)
{
  stopifnot(all(!is.na(spdat)))
  stopifnot(all(coordinates(coords) >= 0 & coordinates(coords) <= 1))

  rl <- list()
  for(i in 1:ncol(spdat))
    rl[[paste0("var", i)]] <- raster::raster(as.matrix(spdat[,i]))
  rl <- raster::stack(rl)

  fn(rl, coords)
}

# experimental, using Matrix
# doesn't seem like this improves speed ... could try it with supernodel?

DEPR_radish_conductance_surface <- function(covariates, coords, directions=8, saveStack=TRUE)
{
  stopifnot(class(covariates) == "RasterStack"  )
  stopifnot(    class(coords) == "SpatialPoints")

  # throw a warning if missing cells are not identical across layers
  spdat   <- as.matrix(raster::getValues(covariates))
  missing <- is.na(rowSums(spdat))
  if (!all(apply(spdat[missing,,drop=FALSE], 1, function(x) all(is.na(x)))))
    warning("Missing cells are not identical across rasters: model selection may be untrustworthy")

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

  # TODO: it might be preferable to form Laplacian matrix Q here?
  # could also do symbolic Cholesky
  out <- list("demes"        = cells,
              "x"            = spdat,
              "adj"          = adj,
              "covariates"   = colnames(spdat), 
              "stack"        = if(saveStack) covariates else NULL)
  class(out) <- "DEPR_radish_graph"
  out
}

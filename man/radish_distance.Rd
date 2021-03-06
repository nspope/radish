% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radish_grid.R
\name{radish_distance}
\alias{radish_distance}
\title{Resistance distances from a parameterized conductance surface}
\usage{
radish_distance(
  theta,
  formula,
  data,
  conductance_model = radish::loglinear_conductance,
  conductance = TRUE,
  covariance = FALSE
)
}
\arguments{
\item{theta}{A matrix of dimension (grid size) x (number of parameters)}

\item{formula}{A formula with the name of a matrix of observed genetic distances on the lhs, and covariates in the creation of \code{data} on the rhs}

\item{data}{An object of class \code{radish_graph} (see \code{\link{conductance_surface}})}

\item{conductance_model}{A function of class \code{radish_conductance_model_factory} (see \code{\link{radish_conductance_model_factory}})}

\item{conductance}{If \code{TRUE}, edge conductance is the sum of cell conductances; otherwise edge conductance is the inverse of the sum of cell resistances (unused; TODO)}

\item{covariance}{If \code{TRUE}, instead of a matrix of resistance distances, return the associated submatrix of the generalized inverse of graph Laplacian}
}
\value{
An object of class \code{radish_grid}
}
\description{
Calculates resistance distances associated with a parameterized conductance surface across
a grid of parameter values.
}
\examples{

library(raster)

data(melip)

covariates <- raster::stack(list(altitude = raster::scale(melip.altitude), 
                                 forestcover = raster::scale(melip.forestcover)))

theta <- as.matrix(expand.grid(forestcover=seq(-1,1,length.out=21), 
                               altitude=seq(-1,1,length.out=21)))

distances <- radish_distance(theta, ~forestcover + altitude, 
                             surface, radish::loglinear_conductance)

ibd <- which(theta[,1] == 0 & theta[,2] == 0)
plot(distances$distance[,,ibd], melip.Fst, pch = 19, 
     xlab = "Null resistance distance (IBD)", ylab = "Fst")

}

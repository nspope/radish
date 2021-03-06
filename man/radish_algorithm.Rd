% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radish_algorithm.R
\name{radish_algorithm}
\alias{radish_algorithm}
\title{Likelihood of parameterized conductance surface}
\usage{
radish_algorithm(
  f,
  g,
  s,
  S,
  theta,
  nu = NULL,
  objective = TRUE,
  gradient = TRUE,
  hessian = TRUE,
  partial = TRUE,
  nonnegative = TRUE,
  validate = FALSE
)
}
\arguments{
\item{f}{A function of class 'conductance_model'}

\item{g}{A function of class 'measurement_model'}

\item{s}{An object of class 'radish_graph'}

\item{S}{A matrix of observed genetic distances}

\item{theta}{Parameters for conductance surface (e.g. inputs to 'f')}

\item{nu}{Number of genetic markers (potentially used by 'g')}

\item{objective}{Compute negative loglikelihood?}

\item{gradient}{Compute gradient of negative loglikelihood wrt theta?}

\item{hessian}{Compute Hessian matrix of negative loglikelihood wrt theta?}

\item{partial}{Compute partial derivatives of negative loglikelihood wrt theta and spatial covariates/observed genetic distances}

\item{nonnegative}{Force regression-like 'measurement_model' to have nonnegative slope?}

\item{validate}{Numerical validation via 'numDeriv' (very slow, use for debugging small examples)}
}
\value{
A list containing at a minimum:
 \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
Additionally, if 'objective == TRUE':
 \item{objective}{(if 'objective') the negative loglikelihood}
 \item{phi}{(if 'objective') fitted values of the nuisance parameters of 'g'}
 \item{boundary}{(if 'objective') is the solution on the boundary (e.g. no genetic structure)?}
 \item{fitted}{(if 'objective') matrix of expected genetic distances among target vertices}
 \item{gradient}{(if 'gradient') gradient of negative loglikelihood with respect to theta}
 \item{hessian}{(if 'hessian') Hessian matrix of the negative loglikelihood with respect to theta}
 \item{partial_X}{(if 'partial') Jacobian of the gradient with respect to the spatial covariates}
 \item{partial_S}{(if 'partial') Jacobian of the gradient with respect to the observed genetic distances}
}
\description{
Calculates likelihood, gradient, hessian, and partial derivatives of a
parameterized conductance surface, given a function mapping spatial data to
conductance and a function mapping resistance distance (covariance) to
genetic distance; using the algorithm in Pope (in prep).
}
\examples{
library(raster)

data(melip)

covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)

radish_algorithm(radish::loglinear_conductance, radish::leastsquares, surface, 
                 ifelse(melip.Fst < 0, 0, melip.Fst), nu = 1000, theta = c(-0.3, 0.3))

}
\references{
Pope NS. In prep. Fast gradient-based optimization of resistance surfaces.
}

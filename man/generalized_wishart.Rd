% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generalized_wishart.R
\name{generalized_wishart}
\alias{generalized_wishart}
\title{Generalized Wishart distance regression}
\usage{
generalized_wishart(
  E,
  S,
  phi,
  nu,
  gradient = TRUE,
  hessian = TRUE,
  partial = TRUE,
  nonnegative = TRUE,
  validate = FALSE
)
}
\arguments{
\item{E}{A submatrix of the generalized inverse of the graph Laplacian (e.g. a covariance matrix)}

\item{S}{A matrix of observed genetic distances}

\item{phi}{Nuisance parameters (see details)}

\item{nu}{Number of genetic markers}

\item{gradient}{Compute gradient of negative loglikelihood with regard to \code{phi}?}

\item{hessian}{Compute Hessian matrix of negative loglikelihood with regard to \code{phi}?}

\item{partial}{Compute second partial derivatives of negative loglikelihood with regard to \code{phi}, \code{E}, \code{S}?}

\item{nonnegative}{Unused}

\item{validate}{Numerical validation via package \code{numDeriv} (very slow, use for debugging small examples)}
}
\value{
A list containing:
 \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
 \item{objective}{(if \code{objective}) the negative loglikelihood}
 \item{fitted}{((if \code{objective}) a matrix of expected genetic distances among target vertices}
 \item{boundary}{(if \code{objective}) is the MLE on the boundary (e.g. no genetic structure)?}
 \item{gradient}{(if \code{gradient}) gradient of negative loglikelihood with respect to phi}
 \item{hessian}{(if \code{hessian}) Hessian matrix of the negative loglikelihood with respect to phi}
 \item{gradient_E}{(if \code{partial}) gradient with respect to the generalized inverse of the graph Laplacian}
 \item{partial_E}{(if \code{partial}) Jacobian of \code{gradient_E} with respect to phi}
 \item{partial_S}{(if \code{partial}) Jacobian of \code{gradient} with respect to S}
 \item{jacobian_E}{(if \code{partial}) a function used for reverse algorithmic differentiation}
 \item{jacobian_S}{(if \code{partial}) a function used for reverse algorithmic differentiation}
}
\description{
A function of class \code{measurement_model} that calculates likelihood,
gradient, hessian, and partial derivatives of nuisance parameters and the
Laplacian generalized inverse, using the generalized Wishart model described in
McCullagh (2009), Peterson et al (2019).
}
\details{
The nuisance parameters are the scaling of the generalized inverse of the graph Laplacian ("tau"; can be zero) and 
a log scalar multiple of the identity matrix that is added to the generalized inverse ("sigma").

TODO: formula
}
\examples{

library(raster)

data(melip)

covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))
surface <- conductance_surface(covariates, melip.coords, directions = 8)

# inverse of graph Laplacian at null model (IBD) 
laplacian_inv <- radish_distance(theta = matrix(0, 1, 2), 
                                 formula = ~forestcover + altitude,
                                 data = surface,
                                 radish::loglinear_conductance, 
                                 covariance = TRUE)$covariance[,,1]

generalized_wishart(laplacian_inv, melip.Fst, nu = 1000, phi = c(0.1, -0.1))

}
\references{
McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica Sinica 19
Peterson et al. 2019. TODO
}
\seealso{
\code{\link{radish_measurement_model}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genetic_distance.R
\name{cov_from_biallelic}
\alias{cov_from_biallelic}
\title{Covariance from allele counts}
\usage{
cov_from_biallelic(Y, N)
}
\arguments{
\item{Y}{A matrix containing allele counts, of dimension (number of demes) x (number of loci)}

\item{N}{A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)}
}
\value{
A matrix containing pairwise covariance
}
\description{
Calculates covariance from counts of the derived allele across biallelic markers,
using normalized allele frequencies
}
\details{
TODO formula
}

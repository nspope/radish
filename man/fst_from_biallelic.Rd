% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genetic_distance.R
\name{fst_from_biallelic}
\alias{fst_from_biallelic}
\title{Fst from allele counts}
\usage{
fst_from_biallelic(Y, N)
}
\arguments{
\item{Y}{A matrix containing allele counts, of dimension (number of demes) x (number of loci)}

\item{N}{A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)}
}
\value{
A matrix containing pairwise Fst
}
\description{
Calculates Fst using the estimator of Bhatia et al (2015)
from counts of the derived allele across biallelic markers
}
\references{
Bhatia et al 2015 TODO
}

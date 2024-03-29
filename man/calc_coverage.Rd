% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{calc_coverage}
\alias{calc_coverage}
\title{Calculate coverage of a given matrix}
\usage{
calc_coverage(cov, min_cov = 1)
}
\arguments{
\item{cov}{integer matrix}

\item{min_cov}{minimum value to count base as covered}
}
\value{
data frame with the following columns:
\itemize{
\item \code{sample} sample names
\item \code{min_cov} minimum coverage threshold (same as \code{min_cov} parameter)
\item \code{covered} number of covered bases
\item \code{total} number of bases
}
}
\description{
This function expects integer matrix where columns are samples and rows are
bases. It normally inputs \code{\link[bsseq:getCoverage]{bsseq::getCoverage()}} result object.
}
\examples{
library(bsseqData)
data(BS.cancer.ex)
getCoverage(BS.cancer.ex) \%>\%
  calc_coverage(min_cov = 10)
}

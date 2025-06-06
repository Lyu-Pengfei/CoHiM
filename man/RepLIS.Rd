% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RepLIS.R
\name{RepLIS}
\alias{RepLIS}
\title{Run RepLIS for rLIS calculation in pairwise replicability analysis}
\usage{
RepLIS(pa, pb, oracle = FALSE, pi, A, f1, f2)
}
\arguments{
\item{pa}{Vector of p-values from study A}

\item{pb}{Vector of p-values from study B}

\item{oracle}{Logical. If TRUE, use oracle inputs; otherwise run EM.}

\item{pi}{Initial vector for pi (oracle mode)}

\item{A}{Initial matrix for A (oracle mode)}

\item{f1}{Vector of f1 values (oracle mode)}

\item{f2}{Vector of f2 values (oracle mode)}
}
\value{
A list with repLIS, fdr, and HMM parameters
}
\description{
Run RepLIS for rLIS calculation in pairwise replicability analysis
}

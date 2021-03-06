% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\docType{data}
\name{dm.fit.eso}
\alias{dm.fit.eso}
\title{Esophageal microbiome data}
\format{A list with five objects}
\description{
A list with 5 objects, including \code{gamma}, the parameter vector for a
Dirichlet distribution. Estimated from the esophageal data set described in
Pei Z et al.: Bacterial biota in the human distal esophagus. Proc Natl Acad Sci U S A 2004,
101:4250-4255. Estimated using package \code{dirmult}.
}
\details{
\itemize{
   \item gamma. Parameter vector for Dirichlet distribution
   \item pi. Mean vector for Dirichlet distribution
   \item loglik. Final log-likelihood value
   \item ite. Number of iterations used
   \item theta. Estimated theta-value
}
}


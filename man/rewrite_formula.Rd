% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergm_estimation_functions.R
\name{rewrite_formula}
\alias{rewrite_formula}
\title{Function to rewrite a formula by replacing the response variable with a specified network
This is mostly for internal use.}
\usage{
rewrite_formula(formula, network)
}
\arguments{
\item{formula}{The original formula to be rewritten. Must be `ergm`-formula.}

\item{network}{The network to be used as the new response variable.}
}
\value{
The rewritten formula with the network as the response variable
}
\description{
Function to rewrite a formula by replacing the response variable with a specified network
This is mostly for internal use.
}

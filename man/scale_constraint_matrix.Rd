% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{scale_constraint_matrix}
\alias{scale_constraint_matrix}
\title{Scale Linear Constraint Matrix}
\usage{
scale_constraint_matrix(L, xs, ys = 1)
}
\arguments{
\item{L}{a matrix giving the linear constraints.}

\item{xs}{a vector of length \code{ncol(L)} giving the scaling
of the model matrix.}

\item{ys}{a double giving the scaling of the response.}
}
\description{
Auxiliary function to scale the linear constraint matrices
to be consistent with the scaled model matrix.
}

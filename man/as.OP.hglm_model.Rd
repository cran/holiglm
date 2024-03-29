% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hglm.R
\name{as.OP.hglm_model}
\alias{as.OP.hglm_model}
\title{Convert to OP}
\usage{
\method{as.OP}{hglm_model}(x)
}
\arguments{
\item{x}{an object inheriting from \code{"hglm_model"}.}
}
\value{
A \pkg{ROI} object of class \code{"OP"}.
}
\description{
Convert an object of class \code{"hglm_model"} into a \pkg{ROI} optimization problem (\code{\link[ROI]{OP}}).
}
\details{
This function is mainly for internal use and advanced users which want of
alter the model object or the underlying optimization problem.
This function converts the model object created by \code{\link{hglm_model}}
into a conic optimization problem solveable via \code{\link[ROI]{ROI_solve}}.
}
\examples{
dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
# Use hglm with option dry_run
model <- hglm(y ~ ., data = dat, dry_run = TRUE)
op <- as.OP(model)
# User hglm_model
x <- model.matrix(y ~ ., data = dat)
model <- hglm_model(x, dat[["y"]])
op <- as.OP(model)
}

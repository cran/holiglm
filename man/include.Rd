% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{include}
\alias{include}
\title{Include Constraint}
\usage{
include(vars)
}
\arguments{
\item{vars}{an integer vector specifying the indices for covariates which have to be in the model.}
}
\value{
A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
}
\description{
Ensures that all covariates specified by \code{vars} coefficients are active.
}
\examples{
dat <- rhglm(100, c(1, 2, 3, 4, 5, 6))
constraints <- include(vars = c("x1", "x3"))
hglm(y ~ ., constraints = constraints, data = dat)
}
\seealso{
Other Constraint-Constructors: 
\code{\link{group_equal}()},
\code{\link{group_inout}()},
\code{\link{group_sparsity}()},
\code{\link{k_max}()},
\code{\link{linear}()},
\code{\link{lower}()},
\code{\link{pairwise_sign_coherence}()},
\code{\link{rho_max}()},
\code{\link{sign_coherence}()},
\code{\link{upper}()}
}
\concept{Constraint-Constructors}

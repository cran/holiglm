% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{rho_max}
\alias{rho_max}
\title{Constraint on the Pairwise Correlation of Covariates}
\usage{
rho_max(
  rho = 0.8,
  exclude = "(Intercept)",
  use = c("everything", "all.obs", "complete.obs", "na.or.complete",
    "pairwise.complete.obs"),
  method = c("pearson", "kendall", "spearman")
)
}
\arguments{
\item{rho}{a value in the range [0,1] specifying, the maximum
allowed collinearity between pairs of covariates.}

\item{exclude}{variables to be excluded form the pairwise
correlation constraints (default is \code{"(Intercept)"}).}

\item{use}{an optional character string giving a method for computing
co-variances in the presence of missing values.
The parameter is passed to \code{\link[stats]{cor}},
therefore see \code{\link[stats]{cor}} for more information.}

\item{method}{a character string indicating which correlation coefficient
is to be computed. See \code{\link[stats]{cor}}
for more information.}
}
\value{
A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
}
\description{
Constraint which ensures that only one covariate out of a pair of covariates with a correlation of at least \code{rho}
will be included in the final model.
}
\examples{
beta <- 1:3
Sigma <- cov_matrix(k = length(beta) - 1L, 1, 2, 0.9)
dat <- rhglm(100, beta, sigma = Sigma)
hglm(y ~ ., constraints = rho_max(0.8), data = dat)
}
\seealso{
Other Constraint-Constructors: 
\code{\link{group_equal}()},
\code{\link{group_inout}()},
\code{\link{group_sparsity}()},
\code{\link{include}()},
\code{\link{k_max}()},
\code{\link{linear}()},
\code{\link{lower}()},
\code{\link{pairwise_sign_coherence}()},
\code{\link{sign_coherence}()},
\code{\link{upper}()}
}
\concept{Constraint-Constructors}

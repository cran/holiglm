% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{predict.hglm}
\alias{predict.hglm}
\title{Predict Method for HGLM Fits}
\usage{
\method{predict}{hglm}(object, newdata = NULL, type = c("link", "response"), ...)
}
\arguments{
\item{object}{a fitted object of class inheriting from "hglm".}

\item{newdata}{an optional data frame containing new observations for which
predictions are to be made. If ommitted, the fitted linear predictors are
used.}

\item{type}{the type of predictions to be made. Possible values are
\code{"link"} (default) or \code{"response"}. If \code{"link"}, the
predictions are in the link scale; if \code{"response"}, the predictions
are transformed to the response scale.}

\item{...}{optional arguments currently ignored.}
}
\value{
A vector of predicted values. If \code{type = "link"}, the predicted
 values are in the link scale; if \code{type = "response"}, the predicted
 values are in the response scale.
}
\description{
Obtains predictions from a fitted holistic generalized linear model object.
}
\examples{
dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
fit <- hglm(y ~ ., constraints = k_max(3), data = dat)
pred <- predict(fit)
pred2 <- predict(fit, newdata=dat)
}

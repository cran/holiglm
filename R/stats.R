
has_intercept <- function(x) UseMethod("has_intercept")

has_intercept.glm <- function(x) {
    isTRUE(attr(terms(x), "intercept") == 1L)
}

has_intercept.hglm_model <- function(x) {
    length(x[["variable_categories"]][["intercept"]]) > 0L
}

has_intercept.matrix <- function(x) {
    nm1 <- colnames(x)[1]
    nm1 == "Intercept" || nm1 == "(Intercept)"
}

# Model Sum of Squares (MSS)
MSS <- function(fit) {
    fitted_values <- fitted.values(fit)
    weights <- weights(fit)
    intercept <- has_intercept(fit)
    if (is.null(weights)) {
        if (intercept) {
            sum((fitted_values - mean(fitted_values))^2)
        } else {
            sum(fitted_values^2)
        }
    } else {
        if (intercept) {
            weighted_mean <- sum(weights * fitted_values / sum(weights))
            sum(weights * (fitted_values - weighted_mean)^2)
        } else {
            sum(weights * fitted_values^2)
        }
    }
}


# Residual Sum of Squares (RSS)
RSS <- function(fit) {
    if (is.null(weights)) {
        sum(residuals(fit)^2)
    } else {
        sum(weights(fit) * residuals(fit)^2)
    }
}


r_squared <- function(fit, adjusted = FALSE) {
    mss <- MSS(fit)
    r2 <- mss / (mss + RSS(fit))
    if (adjusted) {
        if (is.null(fit$qr)) stop("hglm object does not have a proper 'qr' component.")
        n <- NROW(fit$qr$qr) - as.integer(has_intercept(fit))
        r2 <- 1 - (1 - r2) * (n / fit$df.residual)
    }
    r2
}


calculate_dispersion <- function(residuals, df_residual, weights = NULL) {
    if (df_residual <= 0) return(NaN)
    if (any(weights == 0)) {
        warning("observations with zero weight not used for calculating dispersion")
    }
    sum((weights * residuals^2)[weights > 0]) / df_residual
}


#
# The following code is a modified version from the 'src/library/stats/R/glm.R' file.
#
#' @noRd
#' @export
summary.hglm <- function(object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...) {
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion)) {
        est.disp <- if (any(object$family$family == c("poisson", "binomial"))) FALSE else TRUE
        if (est.disp) {
            dispersion <- calculate_dispersion(object$residuals, df.r, object$weights)
        } else {
            dispersion <- 1
        }
    }

    ## calculate scaled and unscaled covariance matrix
    aliased <- is.na(coef(object))
    p <- object$rank
    if (isTRUE(p > 0)) {
        p1 <- seq_len(p)
        if (is.null(Qr <- object$qr)) stop("hglm object does not have a proper 'qr' component.")
        coef.p <- object$coefficients[object$coefficients.selected]
        ## calculate correction for SEs 
        Psi <- object$hglm_model$correction_se
        if (is.null(Psi)) {
            covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        } else {
            u <- Qr$qr[p1, p1, drop = FALSE]
            u[lower.tri(u)] <- 0
            invPsiU <- chol2inv(chol(crossprod(u %*% Psi)))
            covmat.unscaled <- Psi %*% tcrossprod(invPsiU, Psi)
        }
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        ## calculate coef table
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p / s.err

        dn <- c("Estimate", "Std. Error")
        if (!est.disp) { # known dispersion
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
        } else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        } else { # df.r == 0
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        }
        df.f <- NCOL(Qr$qr)
    } else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }

    keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                    "null.deviance", "df.null", "iter", "na.action"), names(object), 0L)
    ans <- c(object[keep],
             list(deviance.resid = residuals(object, type = "deviance"), coefficients = coef.table,
                  aliased = aliased, dispersion = dispersion, df = c(object$rank, df.r, df.f),
                  cov.unscaled = covmat.unscaled, cov.scaled = covmat))

    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled / outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.hglm"
    return(ans)
}


#' @noRd
#' @export
summary.hglm.fit <- summary.hglm


#' @noRd
#' @export
print.summary.hglm <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                               signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE),
                                     c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        df <- if ("df" %in% names(x)) x[["df"]] else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L])) {
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
        } else {
            cat("\nCoefficients:\n")
        }
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
        format(x$dispersion), ")\n\n",
        apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                    format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                    " on",
                    format(unlist(x[c("df.null", "df.residual")])),
                    " degrees of freedom\n"),
              1L, paste, collapse = " "),
        sep = "")
    if (nzchar(mess <- naprint(x$na.action))) {
        cat("  (", mess, ")\n", sep = "")
    }
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
        "\n\n", "Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2L), nsmall = 2L, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}


#' Obtain all Active Coefficients
#'
#' The function returns a logical vector which is \code{TRUE} for all
#' active (i.e., non-zero) coefficients in the fitted model and \code{FALSE} otherwise.
#'
#' @param object an object inheriting from \code{"hglm"} or \code{"hglm.fit"}
#'  from which the active coefficients obtained from.
#' @param ... optional arguments currently ignored.
#' @return a logical vector giving the active coefficients.
#' @aliases acoef
#' @rdname active_coefficients
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' fit <- hglm(y ~ ., constraints = k_max(3), data = dat)
#' active_coefficients(fit)
#' @export
active_coefficients <- function(object, ...) UseMethod("active_coefficients")


#' @noRd
#' @export
active_coefficients.hglm <- function(object, ...) {
    object[["coefficients.selected"]]
}


#' @noRd
#' @export
active_coefficients.hglm.fit <- active_coefficients.hglm


#' @rdname active_coefficients
#' @export
acoef <- active_coefficients

#' Extract Model Coefficients
#'
#' The function returns different types of coefficients for a model.
#'
#' @details
#' The types \code{"scaled"} and \code{"unscaled"} refer to the coefficients
#'  of the scaled/unscaled optimization problem. Type \code{"selected"} refers
#'  to the active coefficients in the model \code{\link{active_coefficients}}.
#'
#' @param object an object inheriting from \code{"hglm"} or \code{"hglm.fit"}
#'  from which the coefficients are obtained from.
#' @param type Default value is \code{"unscaled"}. Allowed values are
#'  \code{"unscaled"}, \code{"scaled"} and \code{"selected"}.
#' @param ... optional arguments currently ignored.
#' @return a vector containing the unscaled, scaled or selected coefficients.
#' @rdname coef.hglm
#' @examples
#' dat <- rhglm(1000, 1:3)
#' fit <- hglm(y ~ ., data = dat)
#' coef(fit)
#' coef(fit, type="scaled")
#' coef(fit, type="selected")
#' @export
coef.hglm <- function(object, type = c("unscaled", "scaled", "selected"), ...) {
    type = match.arg(type)
    if (type == "scaled") {
        no_scaling <- length(object[["hglm_model"]][["scaler"]]) && object[["hglm_model"]][["scaler"]] == "off"
        if (no_scaling) warning("No scaling was used.")
        return(object[["coefficients_scaled"]])
    }
    if (type == "selected") return(object[["coefficients.selected"]])
    object[["coefficients"]]
}

#' @noRd
#' @export
coef.hglm.fit <- coef.hglm


# Extract Response from a Model Frame
#
# Extracts the response from a model frame. A simple wrapper around
# \code{\link[stats]{model.response}} which adds the class information
# \code{"model_response"}.
#
# @param data a model frame.
# @param type a character string giving the storage type of the
#             extracted object. Possible values are \code{"any"},
#             \code{"numeric"} and \code{"double"}.
#
# @seealso
# \code{\link[stats]{model.response}}
#
# @return
# Returns the model response.
model_response <- function(data, type = "any") {
    response <- model.response(data = data, type = type)
    class(response) <- c("model_response", class(response))
    response
}


# Build Model Matrix
model_matrix <- function(object, data = environment(object),
                         contrasts.arg = NULL, xlev = NULL, ...) {
    mm <- model.matrix.default(object = object, data = data,
                               contrasts.arg = contrasts.arg, xlev = xlev, ...)
    class(mm) <- c("model_matrix", class(mm))
    mm
}


#' Predict Method for HGLM Fits
#'
#' Obtains predictions from a fitted holistic generalized linear model object.
#'
#' @param object a fitted object of class inheriting from "hglm".
#' @param newdata an optional data frame containing new observations for which
#'  predictions are to be made. If ommitted, the fitted linear predictors are
#'  used.
#' @param type the type of predictions to be made. Possible values are
#'  \code{"link"} (default) or \code{"response"}. If \code{"link"}, the
#'  predictions are in the link scale; if \code{"response"}, the predictions
#'  are transformed to the response scale.
#' @param ... optional arguments currently ignored.
#' @return A vector of predicted values. If \code{type = "link"}, the predicted
#'  values are in the link scale; if \code{type = "response"}, the predicted
#'  values are in the response scale.
#' @rdname predict.hglm
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' fit <- hglm(y ~ ., constraints = k_max(3), data = dat)
#' pred <- predict(fit)
#' pred2 <- predict(fit, newdata=dat)
#' @export
predict.hglm <- function(object, newdata = NULL, type = c("link", "response"), ...) {
    type <- match.arg(type)

    if (is.null(newdata)) {
        mm <- model_matrix(object)
    } else {
        Terms <- delete.response(terms(object))
        m <- model.frame(Terms, newdata, xlev = object$xlevels)
        mm <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    }

    pred <- drop(mm %*% coef(object))
    if (type == "link") pred else family(object)$linkinv(pred)
}


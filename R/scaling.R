# Vector containing c(center, scale_factor, min_max_distance)
scale_factor_minmax <- function(x) {
    x_min <- min(x)
    x_max <- max(x)
    if (x_min == x_max) { # Should only happen for the intercept.
        if (x_min == 0) {
            return(c(0.0, 0, 0))  # Should give an error in the next step.
        } else {
            return(c(0.0, 1 / x_max, 0))
        }
    }
    min_max_distance <- x_max - x_min
    c(x_min, 1 / (min_max_distance), min_max_distance)
}


scale_factor_standardize <- function(x) {
    x_min <- min(x)
    x_max <- max(x)
    x_mean <- mean(x)
    x_sd <- sd(x)
    if (x_min == x_max) { # Should only happen for the intercept.
        if (x_min == 0) {
            return(c(0.0, 0, 0))  # Should give an error in the next step.
        } else {
            return(c(0.0, 1 / x_sd, 0))
        }
    }
    min_max_distance <- x_max - x_min
    c(x_mean, 1 / x_sd, min_max_distance)
}


scale_model_matrix <- function(x, center = TRUE, method = c("minmax", "standardize")) {
    method <- match.arg(method)
    checkmate::assert_logical(center, len = 1L)
    if (method == "minmax") {
        scale_factor <- scale_factor_minmax
    } else {
        scale_factor <- scale_factor_standardize
    }
    xm <- xs <- double(ncol(x))
    for (i in seq_len(ncol(x))) {
        y <- scale_factor(x[, i])
        if (!center) {
            y[1] <- 0.0
        }
        xm[i] <- y[1]
        xs[i] <- y[2]
        if (all(y == c(0, 1))) {
            next()
        }
        x[, i] <- (x[, i] - y[1]) * y[2]
    }
    attr(x, "xm") <- xm
    attr(x, "xs") <- xs
    x
}


scale_response <- function(x, center = TRUE, method = c("minmax", "standardize")) {
    method <- match.arg(method)
    checkmate::assert_logical(center, len = 1L)
    if (method == "minmax") {
        scale_factor <- scale_factor_minmax
    } else {
        scale_factor <- scale_factor_standardize
    }
    y <- scale_factor(x)
    if (y[3] == 0) {
        stop("repsonse has zero variance")
    }
    if (!center) {
        y[1] <- 0.0
    }
    x <- (x - y[1]) * y[2]
    attr(x, "ym") <- y[1]
    attr(x, "ys") <- y[2]
    x
}


#' Scale Linear Constraint Matrix
#'
#' Auxiliary function to scale the linear constraint matrices
#' to be consistent with the scaled model matrix.
#'
#' @param L a matrix giving the linear constraints.
#' @param xs a vector of length \code{ncol(L)} giving the scaling
#'           of the model matrix.
#' @param ys a double giving the scaling of the response.
#'
#' @export
scale_constraint_matrix <- function(L, xs, ys = 1) {
    UseMethod("scale_constraint_matrix")
}


#' @noRd
#' @export
scale_constraint_matrix.matrix <- function(L, xs, ys = 1) {
    checkmate::assert_numeric(xs, len = ncol(L), any.missing = FALSE)
    if (is.null(ys)) {
        ys <- 1
    }
    L %*% diag(ys / xs)
}


#' @noRd
#' @export
scale_constraint_matrix.simple_triplet_matrix <- function(L, xs, ys = 1) {
    checkmate::assert_numeric(xs, len = ncol(L), any.missing = FALSE)
    if (is.null(ys)) {
        ys <- 1
    }
    scaler <- ys / xs
    m <- match(L[["j"]], seq_along(scaler))
    L[["v"]] <- L[["v"]] * scaler[m]
    L
}

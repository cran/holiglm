#
# Define some short cuts.
#
stm <- simple_triplet_matrix

stzm <- simple_triplet_zero_matrix

stdm <- simple_triplet_diag_matrix


is_power <- function(family) {
    grepl("^", family$link, fixed=TRUE)
}


extract_power <- function(family) {
    assert(check_class(family, "family"), check_true(is_power(family)), combine = "and")
    as.numeric(gsub("mu^", "", family$link, fixed=TRUE))
}


factor_binary <- function(x) {
    assert(check_class(x, "factor"), combine = "and")
    x <- factor(x)
    assert(check_true(length(levels(x)) == 2))
    as.integer(x) - 1L
}


loglike_function <- function(family, approx = "") {
    # envir <- parent.frame(1)
    envir <- parent.env(environment())
    tryCatch({
        if (is_power(family)) {
            mu <- extract_power(family)
            qassert(mu, "R1(0,1)")
            get(sprintf("loglike_%s_power%s", family$family, approx), envir = envir)
        } else {
            get(sprintf("loglike_%s_%s%s", family$family, family$link, approx), envir = envir)
        }
    }, error = function(e) NULL)
}


modify_environment <- function(x, val) {
    assert_character(names(val), any.missing = FALSE)
    for (key in names(val)) {
        x[[key]] <- val[[key]]
    }
}


is_zero_one <- function(x) abs(max(round(x) - x)) < 1e-16


#' Aggregate Binomial Data
#'
#' A simple function for aggregating binomial data, from
#' a form where \code{y} contains only \code{0} and \code{1}
#' and \code{X} could contain duplicated rows, into a format where
#' \code{y} is the matrix of counted successes and failures
#' and \code{X} does not contain duplicates. If \code{X} contains
#' factor variables, the model matrix corresponding to \code{X}
#' will be returned.
#'
#' @param formula a formula object defining the aggregation.
#' @param data a \code{data.frame} to be aggregated.
#' @param as_list a logical giving if the return value
#'          should be a \code{list}. If \code{FALSE} the return value
#'          is a \code{data.frame}.
#' @return A \code{list} (or \code{data.frame}) containing aggregated
#'  binomial data with counted successes and failures.
#' @examples
#' set.seed(12345)
#' data <- data.frame(y = rbinom(50, 1, 0.7),
#'                    a = factor(sample(c(1, 2), 50, TRUE)),
#'                    b = factor(sample(c(1, 2, 3), 50, TRUE)))
#' agg_binomial(y ~ ., data)
#' @export
agg_binomial <- function(formula, data, as_list = TRUE) {
    Y <- model.response(model.frame(data), "any") 
    stopifnot(NCOL(Y) == 1L)
    X <- model.matrix(formula, data)
    i <- do.call(order, as.list(as.data.frame(X)))
    x <- X[i,]
    y <- Y[i]
    b <- duplicated(x)
    group_index <- cumsum(!b)
    x <- x[!b,]
    y <- unclass(table(group_index, y))
    names(dimnames(y)) <- NULL
    rownames(x) <- NULL
    rownames(y) <- NULL
    colnames(y) <- c("failure", "success")
    if (as_list) {
        list(x = x, y = y[, c(2, 1)])    
    } else {
        cbind(as.data.frame(y), as.data.frame(x)[,colnames(x) != "(Intercept)"])
    }
}

pad_TRUE <- function(x, n=0L) c(x, !logical(n))
pad_FALSE <- function(x, n=0L) c(x, logical(n))

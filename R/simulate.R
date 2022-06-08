
test_all_equal <- function(x) UseMethod("test_all_equal")

test_all_equal.integer <- function(x) all(diff(x) == 0L)

check_all_equal <- function(x) UseMethod("check_all_equal")

check_all_equal.integer <- function(x) check_true(test_all_equal(x))

test_all_equal_length <- function(...) {
    test_all_equal(sapply(list(...), length))
}


check_all_equal_length <- function(...) {
    check_true(test_all_equal_length(...))
}


#' Construct Covariance matrix
#'
#' Utility function for constructing covariance matrices based on 
#' a simple triplet format (\code{\link[slam]{simple_triplet_matrix}}).
#' 
#' @param k an integer giving the number of rows and columns of the constructed 
#'          covariance matrix.
#' @param i an integer vector giving the row indices.
#' @param j an integer vector giving the row indices.
#' @param v a numeric vector giving the corresponding values.
#' @return A dense \code{matrix} of covariances.
#' @examples
#' cov_matrix(5, c(1, 2), c(2, 3), c(0.8, 0.9))
#' @export
cov_matrix <- function(k, i, j, v) {
    assert(check_all_equal_length(i, j, v))
    covm <- diag(k)
    for (it in seq_along(v)) {
        covm[j[it], i[it]] <- covm[i[it], j[it]] <- v[it]
    }
    covm
}


#' Random HGLM Data
#'
#' A simple data generator for testing and example purposes.
#'
#' @param n the number of observations to be created.
#' @param beta a numeric vector giving the magnitude of the coefficients
#'      (the first element is assumed to be the intercept).
#' @param sigma a positive-definite symmetric matrix giving the covariance
#'      structure of the covariates (passed to \code{MASS::mvrnorm}).
#' @param family the family of the inverse link.
#' @param truncate_mu a logical giving if mu should be truncated if necessary.
#' @param as_list a logical (default is \code{FALSE}), if \code{TRUE} a \code{list} is returned
#'      otherwise a \code{data.frame} is returned.
#' @param ... addtional optional parameters. The arguments are passed to the
#'  random variables generating function of the response.
#' @return A \code{data.frame} (or \code{list}) containing the generated data.
#' @examples
#' rhglm(10, 1:5)
#' @export
rhglm <- function(n, beta, sigma = diag(length(beta) - 1L), family = gaussian(), truncate_mu = FALSE,
                  as_list = FALSE, ...) {
    assert(check_all_equal(c(length(beta) - 1L, NROW(sigma), NCOL(sigma))))
    assert_integerish(n, lower = 1L, len = 1L, any.missing = FALSE)
    assert_numeric(beta, any.missing = FALSE)
    assert_numeric(sigma, any.missing = FALSE)
    family <- get_family(family)
    family_name <- canonicalize_family_name(family$family)
    k <- length(beta) - 1L
    x <- MASS::mvrnorm(n, mu = double(k), Sigma = sigma)
    # FIXME: If coeficients are provided which in combination with the coefs do not
    #        result in a mu which fullfills the requirements we should change it.
    eta <- as.vector(cbind(1, x) %*% beta)
    mu <- family$linkinv(eta)
    if (isTRUE(truncate_mu)) {
        mu <- bound_mu(mu, family_name)
        truncated_beta <- c(MASS::ginv(cbind(1, x)) %*% family$linkfun(mu))    
    }
    y <- sim_response(family_name, n, mu)
    colnames(x) <- sprintf("x%s", seq_len(NCOL(x)))
    if (as_list) {
        dat <- list(y = y, x = cbind(intercept = 1, x))
    } else {
        dat <- cbind(y = y, as.data.frame(x))
    }
    attr(dat, "mu") <- mu
    attr(dat, "true_beta") <- beta
    if (truncate_mu) {
        attr(dat, "trunc_beta") <- truncated_beta
    }
    dat
}


bound_mu <- function(mu, family_name, eps = .Machine[["double.eps"]]) {
    if (family_name == "binomial") {
        mu[mu >= 1] <- 1 - eps
    } else if (family_name %in% c("Gamma", "poisson", "inverse.gaussian", "negative.binomial")) {
        mu[mu <= 0] <- 0
    }
    mu
}


canonicalize_family_name <- function(family_name) {
    if (startsWith(family_name, "Negative Binomial")) {
        "negative.binomial"
    } else {
        family_name
    }
}


sim_gaussian_response <- function(n, mu, sd = 1, ...) {
    assert_numeric(sd, lower = .Machine[["double.eps"]], len = 1L, any.missing = FALSE)
    rnorm(n = n, mean = mu, sd = sd)
}

sim_binomial_response <- function(n, mu, ...) {
    assert_numeric(mu, any.missing = FALSE, lower = 0, upper = 1)
    rbinom(n, 1, mu)
}

sim_gamma_response <- function(n, mu, scale = 1, ...) {
    assert_numeric(mu, any.missing = FALSE, lower = .Machine[["double.eps"]])
    assert_numeric(scale, any.missing = FALSE, lower = .Machine[["double.eps"]])
    rgamma(n = n, shape = mu / scale, scale = scale)
}

sim_poisson_response <- function(n, mu, ...) {
    assert_numeric(mu, any.missing = FALSE, lower = .Machine[["double.eps"]])
    rpois(n = n, lambda = mu)
}

sim_inverse.gaussian_response <- function(n, mu, lambda = 1, dispersion = 1, ...) {
    assert_numeric(mu, any.missing = FALSE, lower = .Machine[["double.eps"]])
    if (dispersion != 1) {
        lambda <- 1 / dispersion
    }
    SuppDists::rinvGauss(n = n, nu = mu, lambda = lambda)
}

sim_negative.binomial_response <- function(n, mu, theta = 1) {
    assert_numeric(mu, any.missing = FALSE, lower = .Machine[["double.eps"]])
    MASS::rnegbin(n = n, mu = mu, theta = theta)
}


sim_response <- function(family_name, n, mu, ...) {
    fnames <- c("gaussian", "binomial", "Gamma", "poisson")
    assert_choice(family_name, fnames)
    kwargs <- list(...)
    simulate <- list(
        gaussian = sim_gaussian_response,
        binomial = sim_binomial_response,
        Gamma = sim_gamma_response,
        poisson = sim_poisson_response)[[family_name]]
    simulate(n, mu, ...)
}


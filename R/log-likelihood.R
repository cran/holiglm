log1exp <- function(xi, j, n_y_is_0) {
    k <- length(xi)
    si <- c(rep.int(1, k),     3,                4,     6)
    sj <- c(seq_along(xi), k + j, k + n_y_is_0 + j, k + j)
    sv <- c(          -xi,    -1,               -1,     1)
    simple_triplet_matrix(si, sj, sv, nrow = 6L, ncol = (k + 2 * n_y_is_0))
}


check_data <- function(x, y) {
    assert(check_vector(y), check_true(length(y) == NROW(x)), combine = "and")
}


#
# Binomial
#
# TODO:
#   The likelihood should store the variables names, so we
#   later know on which part of the matrices we have to set the constriants.
# @title Binomial Log Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y a numeric vector giving the response variable.
# @param eps a numeric giving the epsilon for the linear constraint.
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_binomial_log <- function(x, y, weights = rep.int(1L, NROW(x)), eps = 1e-7, ...) {
    wy <- y * weights
    assert_integerish(wy, any.missing = FALSE)
    if (is_zero_one(y)) {
        yx <- y %*% x
        y_is_0 <- y == 0L
        weights <- rep.int(1, sum(y_is_0))
    } else {
        yx <- round(wy) %*% x
        y0_count <- round((1 - y) * weights)
        y_is_0 <- y0_count > 0L
        weights <- y0_count[y_is_0]
    }
    n_y_is_0 <- sum(y_is_0)
    op <- OP(-c(beta=yx, double(n_y_is_0), weights), maximum = FALSE)
    L1 <- cbind(x, stzm(nrow(x), 2 * n_y_is_0))
    L2 <- mapply(log1exp, split(x[y_is_0,], seq_len(n_y_is_0)),
                 seq_len(n_y_is_0), MoreArgs = list(n_y_is_0 = n_y_is_0),
                 SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-eps, nrow(x)), rep(rhs, n_y_is_0))
    cones <- c(K_lin(nrow(x)), K_expp(2 * n_y_is_0))
    L <- do.call(rbind, c(list(L1), L2))
    constraints(op) <- C_constraint(L, cones, rhs)
    bounds(op) <- V_bound(ld = -Inf, nobj = length(objective(op)))
    op
}


# @title Binomial Logit Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_binomial_logit <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if (is.null(weights)) {
        yx <- y %*% x
        weights <- rep.int(1L, NROW(x))
    } else {
        yx <- ((weights * y) %*% x)
    }
    m <- nrow(x); n <- ncol(x)
    i <- 3 * seq_len(m) - 2
    # op <- OP(c(beta=-(y %*% x), rep.int(1, m), double(m)), maximum = FALSE)
    op <- OP(c(beta=-yx, weights, double(m)), maximum = FALSE)
    C11 <- stm(rep(i, n), rep(seq_len(n), each = m), -drop(x), 3 * m, n)
    C12 <- stm(i, seq_len(m), rep.int(1, m), 3 * m, m)
    C13 <- stm(i + 2, seq_len(m), rep.int(-1, m), 3 * m, m)
    C1 <- cbind(C11, C12, C13)
    C2 <- cbind(stzm(3 * m, n), C12, -C13)
    C <- rbind(C1, C2)
    cones <- K_expp(2 * m)
    rhs <- c(rep(c(0, 1, 1), m), rep(c(0, 1, 0), m))
    constraints(op) <- C_constraint(C, cones, rhs)
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(C))
    op
}


# @title Binomial Logit Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @references
# K.D. Tocher:The  art  of  simulation. (1963) English Universities Press,London.
# @export
# NOTE: We use the Kocher approximation 
loglike_binomial_probit <- loglike_binomial_logit


#
# Poisson
#
# @title Poisson Log Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_poisson_log <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if (is.null(weights)) {
        yx <- y %*% x
        weights <- rep.int(1L, NROW(x))
    } else {
        yx <- (weights * y) %*% x
    }
    m <- nrow(x); n <- ncol(x)
    i <- 3 * seq_len(m) - 2
    op <- OP(c(-yx, weights))
    A <- cbind(stm(rep(i, n), rep(seq_len(n), each = m), -drop(x), 3 * m, n),
               stm(i + 2, seq_len(m), rep.int(-1, m), 3 * m, m))
    rhs <- rep(c(0, 1, 0), m)
    cones <- K_expp(m)
    constraints(op) <- C_constraint(A, cones, rhs)
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(A))
    op
}


# @title Poisson Identity Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_poisson_identity <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if (is.null(weights)) {
        wx <- x
        wy <- y
    } else {
        wx <- diag(weights) %*% x
        wy <- weights * y
    }
    m <- nrow(x); n <- ncol(x)
    i <- 3 * seq_len(m) - 2
    op <- OP(c(colSums(wx), -wy))
    A <- cbind(stm(rep(i+2, n), rep(seq_len(n), each = m), -drop(x), 3 * m, n),
               stm(i, seq_len(m), rep.int(-1, m), 3 * m, m))
    rhs <- rep(c(0, 1, 0), m)
    cones <- K_expp(m)
    nonneg_const <- C_constraint(cbind(-wx, stzm(m)), K_lin(m), rep(0, m))
    constraints(op) <- c(C_constraint(A, cones, rhs), nonneg_const)
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(A))
    op
}


# @title Poisson Sqrt Log-Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_poisson_sqrt <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if (is.null(weights)) {
        wx <- x
        wy <- y
    } else {
        wx <- diag(sqrt(weights)) %*% x
        wy <- weights * y
    }
    m <- nrow(x); n <- ncol(x)
    #          X          zeta   delta
    op <- OP(c(double(n), 1,     -2 * wy))
    i <- 3 * seq_len(m) - 2
    A <- rbind(stm(c(1, 2), c(n + 1, n + 1), c(-1, 1), ncol = n + 1 + m),
               stm(rep(seq_len(m), n), rep(seq_len(n), each = m), -2 * drop(wx), ncol = n + 1 + m),
               cbind(stm(rep(i + 2, n), rep(seq_len(n), each = m), -drop(x), 3 * m, n),
                     stm(i, seq_len(m) + 1, rep.int(-1, m), 3 * m, m + 1))
              )
    rhs <- c(1, 1, integer(m), rep(c(0, 1, 0), m))
    cones <- c(K_soc(m + 2), K_expp(m))
    constraints(op) <- C_constraint(A, cones, rhs)
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(A))
    op
}


#
# Gaussian data
#
# @title Gaussian Identity Likelihood
# @param x design matrix of dimension \eqn{n \times p}.
# @param y response variable
# @param weights an optional numeric vector (or \code{NULL}) of prior
#        weights to be used in the fitting process. The length
#        of numeric weights should be equal to \code{NROW(x)}.
# @param solver a character string giving the solver name.
# @param ... optional additional variables, currently not used.
# @return An object of class \code{"OP"}.
# @export
loglike_gaussian_identity <- function(x, y, weights = rep.int(1L, NROW(x)), solver = "auto", ...) {
    SOCP_solvers = c("ecos", "clarabel", "scs")
    if (solver %in% SOCP_solvers || solver=="auto" && "ecos" %in% names(ROI::ROI_registered_solvers())) {
        # FIXME: Use a better algorithm to select the solver.
        return(loglike_gaussian_identity_SOCP(x, y, weights = weights, ...))
    }
    if (is.null(weights)) {
        Q <-  2 * t(x) %*% x
        L <- -2 * t(y) %*% x
    } else {
        Q <-  2 * t(x) %*% diag(weights) %*% x
        L <- -2 * t(y) %*% diag(weights) %*% x
    }
    op <- OP(objective = Q_objective(Q = Q, L = L))
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(x))
    op
}


# SOCP formulation
loglike_gaussian_identity_SOCP <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if (!is.null(weights)) {
        x <- diag(sqrt(weights)) %*% x
        y <- diag(sqrt(weights)) %*% y
    }
    m <- NROW(x); n <- NCOL(x)
    op <- OP(c(beta=double(n), SOC=1))
    A <- rbind(
      stm(i = c(1,2), j = c(n+1, n+1), v = c(-1, 1), ncol = n+1),
      cbind(2*x, stzm(m, 1))
    )
    constraints(op) <- C_constraint(A, c(K_soc(m+2)), rhs = c(1, 1, 2*y))
    bounds(op) <- V_bound(ld = -Inf, nobj = n+1)
    op
}


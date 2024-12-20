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

loglike_binomial_log_approx <- function(x, y, weights = rep.int(1L, NROW(x)), eps = 1e-7, ...) {
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(-16.53, -11.88, -8, -7.01, -4.68, -4.23, -3.12, -2.34, -1.92, -1.58, -1.28, -1.03, -0.82, -0.64, -0.49, -0.25, -0.1)
    }
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
    n <- nrow(x); m <- ncol(x)

    f <- function(v) -log(1-exp(v))
    fd <- function(v) 1/(exp(-v)-1)
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(x)[y_is_0,], stdm(1, n_y_is_0)), geq(n_y_is_0), rep(f(vk)-fd(vk)*vk, n_y_is_0))

    # - x * \beta + -log(1-exp(x * \beta))
    op <- OP(c(-(yx), weights), maximum = FALSE)
    C1 = L_constraint(cbind(drop(x), stzm(n, n_y_is_0)), leq(n), double(n))
    C2 = do.call(rbind, lapply(Vx, tangent))
    constraints(op) <- rbind(C1, C2)
    bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
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

loglike_binomial_logit_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(-8.02, -5.88, -4.74, -3.37, -1.75, -0.65, -0.33, -0.04, 0.37, 1.07, 1.69, 2.74, 4.01, 7.08, 8.82)
        #Vx <- c(-6.37, -4.78, -3.79, -3.03, -2.14, -1.47, -0.87, -0.32, 0.2, 0.6, 1.03, 2.02, 3.26, 4.01, 5.4)
    }
    y[y == 0] <- -1
    if (is.null(weights)) {
      yx <- y * x
      weights <- rep.int(1L, NROW(x))
    } else {
      yx <- (weights * y) * x
    }
    n <- nrow(x); m <- ncol(x)
    f <- function(v) log(1+exp(v))
    fd <- function(v) 1/(exp(-v)+1)
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(-yx), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))
    # build_tangent <- function(Vx) {
    #     mat = as.simple_triplet_matrix(-yx)
    #     mat1 = do.call(rbind, lapply(Vx, function(vk) -fd(vk)*(-yx)))
    #     mat2 = stm(seq(n*length(Vx)), rep(seq(n), length(Vx)), rep(1L, length(Vx)*n))
    #     L_constraint(cbind(mat1, mat2), geq(n*length(Vx)), rep(sapply(Vx, function(vk) f(vk)-fd(vk)*vk), each=n))
    # }
    # tangent_cons <- build_tangent(Vx)
    op <- OP(c(beta=double(m), weights), maximum = FALSE)
    tangent_cons <- do.call(rbind, lapply(Vx, tangent))
    # f'(v1)(v-v1)+f(v1)=-v
    v1_cons <- L_constraint(cbind(drop(yx), stdm(1,n)), geq(n), numeric(n))
    # f'(v4)(v-v4)+f(v4)=0
    v4_cons <- L_constraint(cbind(stzm(n, m), stdm(1,n)), geq(n), numeric(n))    
    constraints(op) <- rbind(tangent_cons, v1_cons, v4_cons)
    bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
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


loglike_binomial_probit_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(-7.79, -6.94, -6.01, -5.04, -4.55, -4.04, -3.3, -2.55, -1.94, -1.32, -0.79, -0.24, 0.27, 0.84, 1.62, 7.11, 8.13)
    }

    y[y == 0] <- -1
    if (is.null(weights)) {
        yx <- y * x
        weights <- rep.int(1L, NROW(x))
    } else {
        yx <- (weights * y) * x
    }
    n <- nrow(x); m <- ncol(x)
    f <- function(v) -log(pnorm(v))
    fd <- function(v) -dnorm(v)/pnorm(v)
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(yx), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))

    op <- OP(c(double(m), rep(1, n)), maximum = FALSE)
    constraints(op) = do.call(rbind, lapply(Vx, tangent))
    bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
    op
}


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
    # minimize -y(x*\beta) + exp(x*\beta)
    op <- OP(c(-yx, weights))
    A <- cbind(stm(rep(i, n), rep(seq_len(n), each = m), -drop(x), 3 * m, n),
               stm(i + 2, seq_len(m), rep.int(-1, m), 3 * m, m))
    rhs <- rep(c(0, 1, 0), m)
    cones <- K_expp(m)
    constraints(op) <- C_constraint(A, cones, rhs)
    bounds(op) <- V_bound(ld = -Inf, nobj = ncol(A))
    op
}

loglike_poisson_log_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(1.04, 1.99, 3.34, 4.25, 5.26, 5.64, 5.98, 6.28, 7.24, 7.66, 8.03, 8.48, 8.91, 9.1, 9.29, 9.54, 9.78)
    }
    
    if (is.null(weights)) {
        yx <- y %*% x
        weights <- rep.int(1L, NROW(x))
    } else {
        yx <- (weights * y) %*% x
    }
    n <- nrow(x); m <- ncol(x)
    op <- OP(c(-yx, weights), maximum=FALSE)
    f <- fd <- exp
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(x), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))

    constraints(op) <- do.call(rbind, lapply(Vx, tangent))
    bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
    op
}


# loglike_poisson_log_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
#     if (is.null(weights)) {
#         yx <- y %*% x
#         weights <- rep.int(1L, NROW(x))
#     } else {
#         yx <- (weights * y) %*% x
#     }
#     n <- nrow(x); m <- ncol(x)
#     op <- OP(c(-yx, weights), maximum=FALSE)
#     f <- function(v) exp(v)
#     fd <- function(v) exp(v)
#     tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(x), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))
#     Vx = seq(0, 10, length.out=50)
#     #Vx = c(0, 1, 2, 3, 4, 5)
#     constraints(op) <- do.call(rbind, lapply(Vx, tangent))
#     bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
#     op
# }



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


loglike_poisson_identity_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(1e3, 2.95, 3.7, 4.32, 5.17, 5.51, 5.82, 6.07, 6.3, 6.72, 7.07, 7.4, 7.65, 7.89, 8.31, 8.67, 9)
    }
    if (is.null(weights)) {
        wx <- x
        wy <- y
    } else {
        wx <- diag(weights) %*% x
        wy <- weights * y
    }
    n <- nrow(x); m <- ncol(x)
    #  x * \beta + y * -\log(x * \beta)
    op <- OP(c(colSums(wx), wy))
    f <- function(x) -log(x)
    fd <- function(x) -1/x
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(x), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))
    constraints(op) <- do.call(rbind, lapply(Vx, tangent))
    bounds(op) <- V_bound(li=seq(m), lb=rep.int(-Inf, m), nobj = length(objective(op)))
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


loglike_poisson_sqrt_approx <- function(x, y, weights = rep.int(1L, NROW(x)), ...) {
    # could rewrite this for non SOCP solver like gurobi
    if ("tangent_points" %in% names(attributes(x))) {
        Vx <- attributes(x)$tangent_points
    } else {
        Vx <- c(0.01, 2.95, 3.7, 4.32, 5.17, 5.51, 5.82, 6.07, 6.3, 6.72, 7.07, 7.4, 7.65, 7.89, 8.31, 8.67, 9)
        # approx works but tangent points are garbage :/
    }

    if (is.null(weights)) {
        wx <- x
        wy <- y
    } else {
        wx <- diag(sqrt(weights)) %*% x
        wy <- weights * y
    }
    n <- nrow(x); m <- ncol(x)
    #          X          zeta   delta
    op <- OP(c(double(m), 1, 2 * wy))
    A <- rbind(stm(c(1, 2), c(m + 1, m + 1), c(-1, 1), ncol = m + 1 + n),
               stm(rep(seq_len(n), m), rep(seq_len(m), each = n), -2 * drop(wx), ncol = m + 1 + n))
    rhs <- c(1, 1, integer(n))
    cones <- K_soc(n + 2)
    f <- function(x) -log(x)
    fd <- function(x) -1/x
    tangent <- function(vk) L_constraint(cbind(-fd(vk)*drop(x), stzm(n, 1), stdm(1, n)), geq(n), rep(f(vk)-fd(vk)*vk, n))
    constraints(op) <- rbind(C_constraint(A, cones, rhs), do.call(rbind, lapply(Vx, tangent)))
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

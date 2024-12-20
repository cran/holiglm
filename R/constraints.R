
match_vars <- function(model, vars) {
    mm <- model$x
    mf <- model$frame
    mf_names <- attr(terms(mf), "term.labels")
    is_intercept <- attr(mm, "assign") == 0
    mm_names <- colnames(mm)
    idx_match <- attr(mm, "assign")[!is_intercept]
    intercept_name <- if (any(is_intercept)) "(Intercept)" else NULL
    mf_to_mm <- data.frame(mf_name = c(intercept_name, mf_names[idx_match]), mm_name = mm_names)
    d <- data.frame(vars = vars)
    df <- merge.data.frame(d, mf_to_mm, by.x = "vars", by.y = "mf_name", all.x = TRUE)
    b <- is.na(df[["mm_name"]])
    df[b, "mm_name"] <- df[b, "vars"]
    mf_to_mm[["mf_col_idx"]] <- attr(mm, "assign")
    mf_to_mm[["mm_col_idx"]] <- seq_along(mm_names)
    df <- merge.data.frame(df, mf_to_mm, by = "mm_name", all.x = TRUE)
    if (anyNA(df[["mm_col_idx"]])) {
        msg <- paste("NAME MISMATCH ERROR",
                     "the variables %s can not be found in the model-frame or model-matrix.",
                    "The available names are:\n%s", sep = "\n")
        mismatch_vars <- deparse(df[["vars"]][is.na(df[["mm_col_idx"]])])
        avil_names <- paste(capture.output(mf_to_mm[, c("mf_name", "mm_name")]), collapse = "\n")
        stop(sprintf(msg, mismatch_vars, avil_names))
    }
    df[, c("vars", "mf_name", "mm_name", "mf_col_idx", "mm_col_idx")]
}


match_kvars <- function(model, kvars) {
    df <- match_vars(model, vars = names(kvars))
    df[["value"]] <- kvars[df[["vars"]]]
    df
}


assert_vars <- function(vars) {
    assert(check_integerish(vars, any.missing = FALSE), check_character(vars, any.missing = FALSE))
}


assert_kvars <- function(kvars) {
    assert_numeric(kvars)
    assert_character(names(kvars), any.missing = FALSE)
}


#' Lower Bound
#'
#' Set a lower bound on the coefficients of specific covariates.
#'
#' @param kvars a named vector giving the lower bounds. The names should correspond to the names
#'        of the covariates in the model matrix. 
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @references
#' McDonald, J. W., & Diamond, I. D. (1990).
#' On the Fitting of Generalized Linear Models with Nonnegativity Parameter Constraints.
#' Biometrics, 46 (1): 201–206.
#' \doi{10.2307/2531643}
#'
#' Slawski, M., & Hein, M. (2013).
#' Non-negative least squares for high-dimensional linear models: Consistency and sparse recovery without regularization.
#' Electronic Journal of Statistics, 7: 3004-3056.
#' \doi{10.1214/13-EJS868}
#' @examples
#' set.seed(0)
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' constraints <- lower(c(x2 = 0, x5 = 1))
#' hglm(y ~ ., constraints = constraints, data = dat)
#'
#' # non-negative least squares
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' constraints <- lower(setNames(double(5), paste0("x", 1:5)))
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
lower <- function(kvars) {
    assert_numeric(kvars, any.missing = FALSE)
    x <- list(kvars = kvars)
    structure(x, class = c("lower_bound", "hglmc"))
}


#' Upper Bound
#'
#' Set a upper bound on the coefficient of specific covariates.
#'
#' @param kvars a named vector giving the upper bounds. The names should correspond to the names
#'        of the covariates in the model matrix. 
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @references
#' McDonald, J. W., & Diamond, I. D. (1990).
#' On the Fitting of Generalized Linear Models with Nonnegativity Parameter Constraints.
#' Biometrics, 46 (1): 201–206.
#' \doi{10.2307/2531643}
#'
#' Slawski, M., & Hein, M. (2013).
#' Non-negative least squares for high-dimensional linear models: Consistency and sparse recovery without regularization.
#' Electronic Journal of Statistics, 7: 3004-3056.
#' \doi{10.1214/13-EJS868}
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' constraints <- upper(c(x1 = 0, x4 = 1))
#' hglm(y ~ ., constraints = constraints, data = dat)
#'
#' @export
upper <- function(kvars) {
    assert_numeric(kvars, any.missing = FALSE)
    x <- list(kvars = kvars)
    structure(x, class = c("upper_bound", "hglmc"))
}


#' Linear Constraint
#'
#' @param L a named vector or matrix defining the linear constraints on the coefficients of the covariates.
#' @param dir a character vector giving the direction of the linear constraints.
#' @param rhs a numeric vector giving the right hand side of the linear constraint.
#' @param on_big_m a logical indicating if the constraint should be imposed on the big-M related binary variables.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @references
#' Lawson, C. L., & Hanson, R. J. (1995).
#' Solving least squares problems. Society for Industrial and Applied Mathematics.
#' Society for Industrial and Applied Mathematics.
#' \doi{10.1137/1.9781611971217}
#' @examples
#' # vector constraint
#' beta <- c(1, -2, 3)
#' dat <- rhglm(100, beta)
#' constraints <- c(linear(c(x1 = 2, x2 = 1), "==", 0), rho_max(1))
#' hglm(y ~ ., data = dat, constraints = constraints)
#'
#' # matrix constraint
#' dat <- rhglm(100, c(1, -2, 3, 4, 5, 6, 7))
#' mat <- diag(2)
#' colnames(mat) <- c("x1", "x5")
#' constraints <- c(linear(mat, c("==", "=="), c(-1, 3)), rho_max(1))
#' hglm(y ~ ., data = dat, constraints = constraints)
#'
#' @export
linear <- function(L, dir, rhs, on_big_m = FALSE) {
    assert_numeric(L, any.missing = FALSE)
    assert_numeric(rhs, any.missing = FALSE)
    assert_true(all(dir %in% c(">=", "==", "<=")))
    assert_true(length(dir) == length(rhs))
    x <- list(L = L, dir = dir, rhs = rhs, on_big_m = on_big_m)
    structure(x, class = c("linear", "hglmc"))
}


# Big-M Constraint for the Covariates
#
# The big-M constraint is used to place an integer upper bound $M$ on all the coefficients of the
# covariates. This constraint is a requirement for many of the other
# constraints.
#
# @param m the big-M parameter, an upper bound for the coefficients, needed for the big-M constraint.
# @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
# @family Constraint-Constructors
# @examples
# big_m(10)
# @export
big_m <- function(m = 1000) {
    assert_numeric(m, len = 1, lower = 0)
    if (identical(m, Inf, FALSE, FALSE, FALSE, FALSE)) return(NULL)
    structure(list(big_m = as.double(m)), class = c("big_m", "hglmc"))
}


#' Constraint on the Number of Covariates
#'
#' Constraint on the maximum number of covariates to be used in the model.
#'
#' @param k an positive integer with \eqn{k \leq k_{max}} giving the maximum number of covariates to be used in the model.
#' @note
#' \itemize{
#'   \item If an intercept is used, the upper bound on \eqn{k_{max} + 1} is given by number of columns of the model matrix.
#'   \item If no intercept is used, the upper bound on \eqn{k_{max}} is given by number of columns of the model matrix.    
#' }
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' hglm(y ~ ., constraints = k_max(3), data = dat)
#' @export
k_max <- function(k) {
    if (identical(k, Inf, FALSE, FALSE, FALSE, FALSE)) return(NULL)
    assert_integerish(k, len = 1, lower = 1, any.missing = FALSE)
    structure(list(k_max = as.integer(k)), class = c("k_max", "hglmc"))
}


#' Constraint on the Pairwise Correlation of Covariates
#'
#' Constraint which ensures that only one covariate out of a pair of covariates with a correlation of at least \code{rho}
#' will be included in the final model.
#'
#' @param rho a value in the range [0,1] specifying, the maximum
#'            allowed collinearity between pairs of covariates.
#' @param exclude variables to be excluded form the pairwise
#'            correlation constraints (default is \code{"(Intercept)"}).
#' @param use an optional character string giving a method for computing
#'            co-variances in the presence of missing values.
#'            The parameter is passed to \code{\link[stats]{cor}},
#'            therefore see \code{\link[stats]{cor}} for more information.
#' @param method a character string indicating which correlation coefficient
#'            is to be computed. See \code{\link[stats]{cor}}
#'            for more information.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' beta <- 1:3
#' Sigma <- cov_matrix(k = length(beta) - 1L, 1, 2, 0.9)
#' dat <- rhglm(100, beta, sigma = Sigma)
#' hglm(y ~ ., constraints = rho_max(0.8), data = dat)
#' @export
rho_max <- function(rho = 0.8, exclude = "(Intercept)",
                    use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"),
                    method = c("pearson", "kendall", "spearman")) {
    use <- match.arg(use)
    method <- match.arg(method)
    assert_numeric(rho, len = 1, lower = 0, upper = 1, any.missing = FALSE)
    if (isTRUE(rho == 1)) return(NULL)
    structure(list(rho_max = rho, exclude = exclude, use = use, method = method), class = c("rho_max", "hglmc"))
}


#' Group Sparsity Constraint
#'
#' Constraint which restricts the number of covariates selected from a specific group.
#'
#' @param vars a vector specifying the indices or names of the covariates to which the group constraint 
#'   shall be applied.
#' @param k an integer giving the maximum number of covariates to be included in the model from the specified group.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' dat <- rhglm(100, c(1, 2, 0, 4, 5, 0))
#' constraints <- group_sparsity(c("x1", "x2", "x5"), 1L)
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
group_sparsity <- function(vars, k = 1L) {
    assert_vars(vars)
    assert_integerish(k, len = 1, lower = 1, any.missing = FALSE)
    x <- list(vars = vars, k_max = k)
    structure(x, class = c("group_sparsity", "hglmc"))
}


#' In-Out Constraint
#'
#' Forces coefficients of the covariates in the specified group to be either all zero or all nonzero.
#'
#' @param vars a vector specifying the indices or names of the covariates to which the constraint shall be applied.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' dat <- rhglm(100, c(1, 2, 3, 4, 5, 6))
#' constraints <- group_inout(c("x1", "x2", "x3"))
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
group_inout <- function(vars) {
    assert_vars(vars)
    x <- list(vars = vars)
    structure(x, class = c("group_inout", "hglmc"))
}


#' Group Equal Constraint
#'
#' Forces all covariates in the specified group to have the same coefficient.
#'
#' @param vars a vector specifying the indices or names of the covariates to which the constraint shall be applied.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' dat <- rhglm(100, c(1, 2, 3, 4, 5, 6))
#' constraints <- group_equal(vars = c("x1", "x3"))
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
group_equal <- function(vars) {
    assert_vars(vars)
    x <- list(vars = vars)
    structure(x, class = c("group_equal", "hglmc"))
}


#' Include Constraint
#'
#' Ensures that all covariates specified by \code{vars} coefficients are active.
#'
#' @param vars an integer vector specifying the indices for covariates which have to be in the model.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @examples
#' dat <- rhglm(100, c(1, 2, 3, 4, 5, 6))
#' constraints <- include(vars = c("x1", "x3"))
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
include <- function(vars) {
    assert_vars(vars)
    x <- list(vars = vars)
    structure(x, class = c("include", "hglmc"))
}


#' Sign Coherence Constraint
#'
#' Constraint which ensures that the coefficients of the specified covariates have
#' a coherent sign.
#'
#' @param vars a character vector giving the names of the covariates the
#'   constraint should be applied to.
#' @param big_m a double giving the big-M parameter.
#' @param eps a double giving the epsilon used to ensure that the constraint holds.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @family Constraint-Constructors
#' @references
#' Carrizosa, E., Olivares-Nadal, A. V., & Ramírez-Cobo, P. (2020).
#' Integer Constraints for Enhancing Interpretability in Linear Regression.
#' SORT. Statistics and Operations Research Transactions, 44: 67-98.
#' \doi{10.2436/20.8080.02.95}.
#' @examples
#' dat <- rhglm(100, c(1, -2, 3, 4, 5, 6))
#' constraints <- sign_coherence(c("x1", "x3"))
#' hglm(y ~ ., constraints = constraints, data = dat)
#' @export
sign_coherence <- function(vars, big_m = 100, eps = 1e-6) {
    x <- list(vars = vars, big_m = big_m, eps = eps)
    structure(x, class = c("sign_coherence", "hglmc"))
}


#' Pairwise Sign Coherence
#'
#' Ensures that coefficients of covariates which exhibit strong pairwise correlation
#' have a coherent sign.
#'
#' @param rho a value in the range [0,1] specifying the maximum
#'            allowed collinearity between pairs of covariates.
#' @param exclude a character vector giving the names of the covariates to be
#'                excluded from the constraint (default is "(Intercept)").
#' @param big_m a double giving the big-M parameter.
#' @param eps a double giving the epsilon for the equal sign constraint.
#'            Since most numerical solvers can only handle constraints up to some epsilon,
#'            e.g., the constraint \eqn{A x \geq b} is typically transformed to
#'            \eqn{|A x - b| \geq 0}. By providing an \code{eps}\eqn{ > 0} and changing
#'            the constraint to \eqn{|A x - b| \geq} \code{eps} we can ensure \eqn{|A x - b| > 0}.
#' @param use an optional character string giving a method for computing
#'            covariances in the presence of missing values.
#'            The parameter is passed to \code{\link[stats]{cor}},
#'            therefore see \code{\link[stats]{cor}} for more information.
#' @param method a character string indicating which correlation coefficient
#'            is to be computed.
#'            The parameter is passed to \code{\link[stats]{cor}},
#'            therefore see \code{\link[stats]{cor}} for more information.
#' @return A holistic generalized model constraint, object inheriting from class \code{"hglmc"}.
#' @examples
#' constraints <- c(k_max(7), pairwise_sign_coherence())
#' @family Constraint-Constructors
#' @references
#' Carrizosa, E., Olivares-Nadal, A. V., & Ramírez-Cobo, P. (2020).
#' Integer Constraints for Enhancing Interpretability in Linear Regression.
#' SORT. Statistics and Operations Research Transactions, 44: 67-98.
#' \doi{10.2436/20.8080.02.95}.
#' @export
pairwise_sign_coherence <- function(
    rho = 0.8,
    exclude = "(Intercept)",
    big_m = 100,
    eps = 1e-6,
    use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"),
    method = c("pearson", "kendall", "spearman")) {
    if (isTRUE(rho == 1)) return(NULL)
    use <- match.arg(use)
    method <- match.arg(method)
    assert_numeric(rho, len = 1, lower = 0, upper = 1, any.missing = FALSE)
    assert_numeric(big_m, len = 1, lower = 0, any.missing = FALSE)
    assert_numeric(eps, len = 1, lower = 0, any.missing = FALSE)
    x <- list(rho_max = rho, exclude = exclude, big_m = big_m, eps = eps, use = use, method = method)
    structure(x, class = c("pairwise_sign_coherence", "hglmc"))
}


#' Generic Functions for \code{hglmc} Objects
#'
#' Generic functions for holistic 'GLM' constraints.
#'
#' The 'HGLM' constraints are all of class \code{"hglmc"}
#' and can be combined with the typical combine function \code{c()}.
#' To verify that an object is a 'HGLM' constraint, the function
#' \code{is.hglmc} can be used.
#'
#' @param ... multiple objects inheriting from \code{"hglmc"} to be combined.
#' @return
#' The combine function \code{c()} returns an object of class \code{"hglmc"}.
#' The \code{is.hglmc} function returns \code{TRUE} if the object inherits
#' from class \code{"hglmc"} otherwise \code{FALSE}.
#' @examples
#' constraints <- c(k_max(7), pairwise_sign_coherence())
#' is.hglmc(constraints)
#' @name hglmc
#' @rdname hglmc
NULL

#' @rdname hglmc
#' @export
c.hglmc <- function(...) {
    x <- list(...)
    assert_true(all(sapply(x, is.hglmc)))
    structure(x, class = "lohglmc") ## list of hglmc
}


#' @param x an R object.
#' @rdname hglmc
#' @export
is.hglmc <- function(x) inherits(x, "hglmc")


#' @noRd
#' @export
print.hglmc <- function(x, ...) {
    writeLines(sprintf("Holistic GLM Sparsity Constraint of Type '%s'", class(x)[1L]))
    str(unclass(x))
}


#' @noRd
#' @export
print.lohglmc <- function(x, ...) {
    writeLines("List of Holistic GLM Sparsity Constraints")
    print(unclass(x))
}



# @title HGLM Constraints
#
# @description Helper function to setup all the desired constraints.
#
# @param ... additional constraints of class "hglmc".
# @param k_max k_max an integer giving the maximum number of covariates to be used.
# @param rho_max a value in the range [0,1] specifying, the maximum
#   allowed collinearity between pairs of covariates.
# @return An object of class \code{"lohglmc"} (list of holistic generalized linear model constraints).
# @rdname hglm_constraints
# @examples
# as.lohglmc(list(k_max(3), rho_max(0.8), group_equal(c("x1", "x3"))))
# c(k_max(3), rho_max(0.8), group_equal(c("x1", "x3")))
hglm_constraints <- function(..., k_max = Inf, rho_max = 0.8) {
    cstrts <- list(...)
    assert_true(all(sapply(cstrts, is.hglmc)))
    cstrts <- c(list(k_max(k_max), rho_max(rho_max)), cstrts)
    structure(cstrts[lengths(cstrts) > 0], class = "lohglmc") ## list of hglmc
}


# @rdname hglm_constraints
# @param x a list of \code{"hglmc"} objects to be coerced into an object of class
#          \code{"lohglmc"} (list of holistic generalized linear model constraints).
# @examples
# as.lohglmc(NULL)
as.lohglmc <- function(x) {
    if (length(x) == 0L) return(structure(list(), class = "lohglmc"))
    if (inherits(x, "lohglmc")) return(x)
    if (inherits(x, "hglmc")) return(structure(list(x), class = "lohglmc"))
    assert_list(x)
    assert_true(all(sapply(x, is.hglmc)))
    structure(x[lengths(x) > 0], class = "lohglmc") ## list of hglmc
}


# Add Constraint to HGLM Model
#
# Add a constraint to a HGLM model (\code{"hglm_model"}).
#
# @param model a HGLM model (object of class \code{"hglm_model"}).
# @param constraint a HGLM constraint constructor (object of class \code{"hglmc"}).
# @family constraint
# @rdname add_constraint
# @return An object of class \code{"hglm_model"}.
# @examples
# beta <- c(1, 2, 0, 4, 5, 0)
# dat <- rhglm(100, beta)
# model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
# model <- add_constraint(model, group_sparsity(c("x1", "x2", "x5"), 1L))
# coef(hglm_fit(model))
# @export
add_constraint <- function(model, constraint) {
    UseMethod("add_constraint", constraint)
}


# @noRd
# @export
add_constraint.big_m <- function(model, constraint) {
    add_bigm_constraint(model, big_m = constraint[["big_m"]])
}


# @noRd
# @export
add_constraint.k_max <- function(model, constraint) {
    add_global_sparsity_constraint(model, k_max = constraint[["k_max"]])
}


# @noRd
# @export
add_constraint.rho_max <- function(model, constraint) {
    add_pairwise_multicollinearity_constraint(model, constraint[["rho_max"]],
        constraint[["exclude"]], constraint[["use"]], constraint[["method"]])
}


# @noRd
# @export
add_constraint.group_sparsity <- function(model, constraint) {
    add_group_sparsity_constraint(model, constraint[["vars"]], constraint[["k_max"]])
}


# @noRd
# @export
add_constraint.group_inout <- function(model, constraint) {
    add_group_inout_constraint(model, constraint[["vars"]])
}


# @noRd
# @export
add_constraint.group_equal <- function(model, constraint) {
    add_group_equal_constraint(model, constraint[["vars"]])
}


# @noRd
# @export
add_constraint.include <- function(model, constraint) {
    add_in_constraint(model, constraint[["vars"]])
}


# @noRd
# @export
add_constraint.lower_bound <- function(model, constraint) {
    add_lower_bound(model, constraint[["kvars"]])
}


# @noRd
# @export
add_constraint.upper_bound <- function(model, constraint) {
    add_upper_bound(model, constraint[["kvars"]])
}


# @noRd
# @export
add_constraint.linear <- function(model, constraint) {
    add_linear_constraint(model, constraint[["L"]], constraint[["dir"]],
                          constraint[["rhs"]], constraint[["on_big_m"]])
}


# @noRd
# @export
add_constraint.sign_coherence <- function(model, constraint) {
    add_sign_constraint(model, constraint[["vars"]], constraint[["big_m"]], constraint[["eps"]])
}


# @noRd
# @export
add_constraint.pairwise_sign_coherence <- function(model, constraint) {
    add_pairwise_multicollinearity_equal_sign_constraint(model, 
        constraint[["rho_max"]], constraint[["exclude"]], constraint[["big_m"]],
        constraint[["eps"]], constraint[["use"]], constraint[["method"]])
}


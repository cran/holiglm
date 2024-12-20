# This function mimics the first part of the stats::glm function.
get_family <- function(family) {
    if (inherits(family, "family")) return(family)
    if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) return(family())
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    family
}


# This function mimics the data handling of glm.
model_data <- function(formula, data, weights) {
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    ## need stats:: for non-standard evaluation
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (!is.null(weights)) {
        ## avoid any problems with 1D or nx1 arrays by as.vector.
        mf[["(weights)"]] <- weights <- as.vector(weights)
    }

    mt <- attr(mf, "terms") # allow model.frame to have updated it

    Y <- model.response(mf, "any") # e.g. factors are allowed
    ## avoid problems with 1D arrays, but keep names
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
    }
    ## null model support
    contrasts <- NULL
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0L)
    list(X = X, Y = Y, weights = weights, mt = mt, mf = mf)
}


#
scale_model_data <- function(md, scale=c("auto", "center_standardization", "center_minmax", "standardization", "minmax", "off"),
    scale_response = FALSE, family="gaussian") {
    # scale x and store the data necessary for the backtransformation in md
    scale <- match.arg(scale)
    intercept <- has_intercept(md$X)
    center <- substr(scale,1,6)=="center"
    if (scale=="standardization" || scale=="center_standardization") {
        if (intercept && center) {
            xm <- c(0, apply(md$X[,-1, drop=FALSE], 2, mean))
            xs <- c(1, apply(md$X[,-1, drop=FALSE], 2, sd))
        } else if (intercept) {
            # no center, but intercept
            xm <- rep(0, NCOL(md$X))
            xs <- c(1, apply(md$X[,-1, drop=FALSE], 2, sd))
        } else{
            if (center) {
                warning("Intercept is deactivated. Data is not shifted towards min. Use scaler='standardization' instead.")
            }
            xm <- rep(0, NCOL(md$X))
            xs <- apply(md$X, 2, sd)
        }
    } else if (scale=="minmax" || scale=="center_minmax" || scale=="auto") {
        if (intercept && center || scale=="auto") {
            xm <- c(0, apply(md$X[,-1, drop=FALSE], 2, min))
            xs <- c(1, apply(md$X[,-1, drop=FALSE], 2, max)) - xm
        } else if (intercept) {
            # no center, but intercept
            xm <- rep(0, NCOL(md$X))
            xs <- c(1, (apply(md$X[,-1, drop=FALSE], 2, max) - apply(md$X[,-1, drop=FALSE], 2, min)))
        } else {
            if (center) {
                warning("Intercept is deactivated. Data is not shifted towards min. Use scaler='minmax' instead.")
            }
            xm <- rep(0, NCOL(md$X))
            xs <- apply(md$X, 2, max) - apply(md$X, 2, min)
        }
    } else {
        xm <- rep(0, NCOL(md$X))
        xs <- rep(1, NCOL(md$X))
    }

    if (scale!="off") {
        md$X <- sweep(md$X, 2L, xm, FUN = "-", check.margin = FALSE)
        md$X <- sweep(md$X, 2L, xs, FUN = "/", check.margin = FALSE)
    }
    attr(md$X, "xm") <- xm
    attr(md$X, "xs") <- xs

    if (scale_response) {
        if (family == "gaussian") {
            ym <- if(intercept) mean(md$Y) else 0
            ys <- sd(md$Y)
            if (abs(ys) <= 1e-4) ys <- 1
        } else if (family == "poisson") { # cater for the extension with additional families
            ym <- 0
            ydiff <- max(md$Y)
            ys <- if(ydiff > 0L && intercept) ydiff else 1
        }
        md$Y <- md$Y - ym
        md$Y <- md$Y / ys
    } else {
        ym <- 0L
        ys <- 1L
    }
    attr(md$Y, "ym") <- ym
    attr(md$Y, "ys") <- ys
    md
}

unscale_model_data <- function(x) {
    xm <- attr(x, "xm")
    xs <- attr(x, "xs")

    if (is.null(xm)) {
        xm <- rep(0, NCOL(x))
    }
    if (is.null(xs)) {
        xs <- rep(1, NCOL(x))
    }

    x <- sweep(x, 2L, xs, FUN = "*", check.margin = FALSE)
    x <- sweep(x, 2L, xm, FUN = "+", check.margin = FALSE)
    x
}

unscale_response <- function(y) {
    ym <- attr(y, "ym")
    ys <- attr(y, "ys")

    if (is.null(ym)) {
        ym <- 0
    }
    if (is.null(ys)) {
        ys <- 1
    }

    y <- (y * ys) + ym
    y
}


get_scale <- function(model) {
    x_scales <- attr(model[["x"]], "xs")
    if (is.null(x_scales)) {
        x_scales <- rep.int(1, NCOL(model[["x"]]))
    }
    x_scales
}

get_scale_response <- function(model) {
    ys <- attr(model[["y"]], "ys")
    if (is.null(ys)) {
        ys <- 1
    }
    ys
}

get_center <- function(model) {
    x_centers <- attr(model[["x"]], "xm")
    if (is.null(x_centers)) {
        x_centers <- rep.int(0, NCOL(model[["x"]]))
    }
    x_centers
}

get_center_response <- function(model) {
    ym <- attr(model[["y"]], "ym")
    if (is.null(ym)) {
        ym <- 0
    }
    ym
}


is_probit <- function(family, approx = FALSE) {
    isTRUE(family$family == "binomial") && isTRUE(family$link == "probit") && !approx
}


approx_inverse <- function(x, family, approx = FALSE) {
    if (is_probit(family, approx)) {
        # K.D. Tocher. 1963. The art of simulation. English Universities Press, London.
        return(x / (2 * sqrt(2 / pi)))
    }
    x
}


# choose_solver("auto", binomial())
# solver <- choose_solver("auto", gaussian())
# choose_solver("auto", gaussian(), FALSE)
# ROI_require_solver()
choose_solver <- function(solver, family, is_mip = NULL, is_conic = FALSE) {
    if (isTRUE(solver == "auto")) {
        if (!is_conic && isTRUE(family[["family"]] == "gaussian") && isTRUE(family[["link"]] == "identity")) {
            solvers <- c("ecos", "gurobi", "cplex", "mosek")
            if (!isTRUE(is_mip) && !is.null(is_mip)) {
                solvers <- c(solvers, "qpoases", "quadprog", "scs", "osqp")
            }
        } else {
            solvers <- c("ecos", "mosek")
        }
        inst_solvers <- names(ROI::ROI_installed_solvers())
        # Since 'ecos' is in depends there is always at least one solver available.
        solver <- solvers[solvers %in% inst_solvers][1L]
    }
    solver
}


#' @title Fitting Holistic Generalized Linear Models
#'
#' @description Fit a generalized linear model under holistic constraints.
#'
#' @details In the case of binding linear constraints the standard errors are corrected,
#'          more information about the correction can be found in
#'          Schwendinger, Schwendinger and Vana (2024) \doi{10.18637/jss.v108.i07}.
#'
#' @param formula an object of class \code{"formula"} giving the symbolic description
#'   of the model to be fitted.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param data a \code{data.frame} or \code{matrix} giving the data for the estimation.
#' @param constraints a list of 'HGLM' constraints stored in a list of class \code{"lohglmc"}.
#'                    Use \code{NULL} to turn off constraints.
#' @param weights an optional vector of 'prior weights' to be used for the estimation.
#' @param scaler a character string giving the name of the scaling function (default is \code{"auto"})
#'               to be employed for the covariates. This typically does not need to be changed.
#' @param scale_response a boolean whether the response shall be standardized or not. Can only
#'                       be used with family \code{gaussian()}. Default is \code{TRUE} for
#'                       family \code{gaussian()} and \code{FALSE} for other families.
#' @param big_m an upper bound for the coefficients, needed for the big-M constraint.
#'  Required to inherit from \code{"hglmc"}. Currently constraints created by
#'  \code{group_sparsity()}, \code{group_inout()},
#'  \code{include()} and \code{group_equal()} use the big-M value specified here.
#' @param solver a character string giving the name of the solver to be used for the estimation.
#' @param control a list of control parameters passed to \code{ROI_solve}.
#' @param dry_run a logical; if \code{TRUE} the model is not fit but only constructed.
#' @param approx a logical; if \code{TRUE} uses linear approximation of log-likelihood.
#' @param object_size a character string giving the object size, allowed values
#'  are \code{"normal"} and \code{"big"}. If \code{"big"} is choosen, also
#'  the \pkg{ROI} solution and the \code{"hglm_model"} object are returned.
#' @param ... For ‘approx’: further arguments passed to or from other methods.
#' @return An object of class \code{"hglm"} inheriting from \code{"glm"}.
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' hglm(y ~ ., constraints = NULL, data = dat)
#' # estimation without constraints
#' hglm(y ~ ., constraints = NULL, data = dat)
#' # estimation with an upper bound on the number of coefficients to be selected
#' hglm(y ~ ., constraints = k_max(3), data = dat)
#' # estimation without intercept
#' hglm(y ~ . - 1, data = dat)
#'
#' @references
#' Schwendinger B., Schwendinger F., Vana L. (2024).
#' Holistic Generalized Linear Models
#' \doi{10.18637/jss.v108.i07}
#'
#' Bertsimas, D., & King, A. (2016).
#' OR Forum-An Algorithmic Approach to Linear Regression
#' Operations Research 64(1):2-16.
#' \doi{10.1287/opre.2015.1436}
#'
#' McCullagh, P., & Nelder, J. A. (2019).
#' Generalized Linear Models (2nd ed.)
#' Routledge.
#' \doi{10.1201/9780203753736}.
#'
#' Dobson, A. J., & Barnett, A. G. (2018).
#' An Introduction to Generalized Linear Models (4th ed.)
#' Chapman and Hall/CRC.
#' \doi{10.1201/9781315182780}
#'
#' Chares, Robert. (2009).
#' “Cones and Interior-Point Algorithms for Structured Convex Optimization involving Powers and Exponentials.”
#'
#' @rdname hglm
#' @name hglm
#' @export
hglm <- function(formula, family = gaussian(), data, constraints = NULL,
                 weights = NULL, scaler = c("auto", "center_standardization",
                 "center_minmax", "standardization", "minmax", "off"),
                 scale_response = NULL, big_m = 100, solver = "auto", control = list(),
                 dry_run = FALSE, approx = FALSE, object_size = c("normal", "big"), ...) {
    tstart <- Sys.time()
    dots <- list(...)
    object_size <- match.arg(object_size)
    constraints <- as.lohglmc(constraints)
    if (length(constraints) == 0L) constraints <- NULL
    if (missing(data)) {
        data <- environment(formula)
    }
    cal <- match.call()
    scaler <- match.arg(scaler)
    family <- get_family(family)
    solver <- choose_solver(solver, family = family, is_mip = TRUE)
    if (is.null(scale_response)) {
        if (family$family == "gaussian") {
            scale_response <- TRUE
        } else {
            scale_response <- FALSE
        }
    }
    if (scale_response && !any(family$family %in% c("gaussian", "poisson"))) {
        msg <- c("Scaling the response is only available for families gaussian() and poisson().\n  ",
                 sprintf("Deactivated scaling response for family %s(%s).", family$family, family$link))
        warning(msg)
        scale_response <- FALSE
    }
    md <- model_data(formula, data, weights)
    if ((attr(md[["mt"]], "intercept") == 0L) && scaler == "auto") {
        scaler <- "off"
    }
    if (is_probit(family, approx)) {
        md$X <- md$X * (2 * sqrt(2 / pi))  # K.D. Tocher. 1963. The art of simulation. English Universities Press, London.
    }

    if ("tangent_points" %in% names(dots)) {
        attr(md$X, "tangent_points") <- dots$tangent_points
    }
    md <- scale_model_data(md, scaler, scale_response = scale_response, family = family$family)
    tstart_model <- Sys.time()
    model <- hglm_model(x = md$X, y = md$Y, family = family, weights = md$weights,
                        frame = md$mf, solver = solver, approx = approx)
    model_build_time <- Sys.time() - tstart_model  # contains building the likelihood
    model[["scaler"]] <- scaler
    fit <- hglm_fit(model, constraints = constraints, big_m = big_m, solver = solver,
                    control = control, dry_run = dry_run, approx = approx, object_size = object_size)
    attr(fit, "object_size") <- object_size
    if (isTRUE(dry_run)) return(fit)
    total_time <- Sys.time() - tstart
    fit$timing <- c(fit$timing, model_build_time = model_build_time, total_time = total_time)
    solver <- fit[["solution_status"]][["msg"]][["solver"]]
    structure(c(fit, list(call = cal, formula = formula, terms = md$mt,
                data = data, constraints = constraints, control = control,
                solver = solver, contrasts = attr(md$X, "contrasts"),
                xlevels = .getXlevels(md$mt, md$mf))),
                class = c("hglm", "glm", "lm"))
}


#' @rdname hglm
#' @name holiglm
#' @export
holiglm <- hglm


#' @param k_seq an integer vector giving the values of \code{k_max} for which the model
#'              should be estimated.
#' @param parallel whether estimation of sequence shall be parallelized
#' @rdname hglm
#' @name hglm_seq
#' @references
#' Chen, J., & Chen, Z. (2008).
#' Extended Bayesian information criteria for model selection with large model spaces.
#' Biometrika, 95 (3): 759–771.
#' Oxford University Press.
#' \doi{10.1093/biomet/asn034}
#'
#' Zhu, J., Wen, C., Zhu, J., Zhang, H., & Wang, X. (2020).
#' A polynomial algorithm for best-subset selection problem.
#' Proceedings of the National Academy of Sciences, 117 (52): 33117–33123.
#' \doi{10.1073/pnas.2014241117}
#'
#' @export
hglm_seq <- function(k_seq, formula, family = gaussian(), data, constraints = NULL,
                     weights = NULL, scaler = c("auto", "center_standardization",
                     "center_minmax", "standardization", "minmax", "off"),
                     big_m = 100, solver = "auto", control = list(), object_size = c("normal", "big"),
                     parallel = FALSE) {
    object_size <- match.arg(object_size)
    scaler <- match.arg(scaler)
    constraints <- as.lohglmc(constraints)
    moma <- model.matrix(formula, data = data)
    p <- NCOL(moma)
    io <- as.integer(any(colnames(moma) == "(Intercept)"))
    if (missing(k_seq)) {
        k_seq <- seq_len(p - io)
    }
    if (max(k_seq) > (p - io)) {
        stop("the model matrix has p = %i columns, therefore max(k_seq) <= (p - 1L) is required", p)
    }
    if (length(constraints) == 0L) {
        constraints <- c(k_max(1))
    }
    if (!"k_max" %in% connames(constraints)) {
        constraints[[length(constraints) + 1L]] <- k_max(1)
    }
    i <- which("k_max" == connames(constraints))
    if (length(i) > 1L) stop("multiple k_max defined")
    fits <- vector("list", length(k_seq))

    stopf = function(fmt, ..., domain = "R-holiglm") stop(gettextf(fmt, ..., domain = domain), domain = NA, call. = FALSE)

    if (parallel) {
        if (!requireNamespace("parallel", quietly = TRUE))
            stopf("Using hglm_seq with parallel=TRUE requires 'parallel' package. Please install 'parallel' or switch off parallel computing.")
        for (j in seq_along(k_seq)) {
            constraints[[i]] <- k_max(k_seq[j])
            fits[[j]] <- parallel::mcparallel({
                fit <- hglm(formula = formula, family = family, data = data, constraints = constraints,
                    weights = weights, scaler = scaler, big_m = big_m, solver = solver,
                    control = control, object_size = object_size)
                fit[["k_max"]] <- k_seq[j]
                fit
            })
        }
    } else {
        for (j in seq_along(k_seq)) {
            constraints[[i]] <- k_max(k_seq[j])
            fit <- hglm(formula = formula, family = family, data = data, constraints = constraints,
                    weights = weights, scaler = scaler, big_m = big_m, solver = solver,
                    control = control, object_size = object_size)
            fit[["k_max"]] <- k_seq[j]
            fits[[j]] <- fit
        }
    }
    if (parallel) {
        fork.res = tryCatch(parallel::mccollect(fits))
        ## check for any errors in FUN, warnings are silently ignored
        fork.err = vapply(fork.res, inherits, FALSE, "try-error", USE.NAMES=FALSE)
        if (any(fork.err))
            stopf("hglm_seq received an error(s) when evaluating FUN:\n%s",
                paste(unique(vapply(fork.res[fork.err], function(err) attr(err,"condition",TRUE)[["message"]], "", USE.NAMES=FALSE)), collapse="\n"))
        fits = unname(fork.res)
    }
    class(fits) <- "hglm_seq"
    fits
}


get_k_max <- function(x) {
    i <- which("k_max" == connames(x$constraints))
    x$constraints[[i]][["k_max"]]
}

# BIC <- function(object, ...) {
#     lls <- logLik(object)
#     -2 * as.numeric(lls) + attr(lls, "df") * log(nobs(lls))
# }

EBIC <- function(object, ..., g = 1) {
    assert_numeric(g, lower=0, upper=1)
    lls <- logLik(object)
    -2 * as.numeric(lls) + log(nobs(lls)) * attr(lls, "df") + 2 * g * log(lchoose(length(coef(object)), attr(lls, "df")-2L))
}

SIC <- function(object, ...) {
    lls <- logLik(object)
    n <- nobs(object)
    n * as.numeric(lls) + log(NCOL(attr(object$terms, "factors"))) * log(log(n)) * (attr(lls, "df")-1L)
}


#' @noRd
#' @export
print.hglm_seq <- function(x, ...) {
    k_max <- sapply(x, get_k_max)
    aic <- round(sapply(x, AIC), 2L)
    bic <- round(sapply(x, BIC), 2L)
    ebic <- round(sapply(x, EBIC), 2L)
    sic <- round(sapply(x, SIC))
    d <- cbind(k_max = k_max, aic = aic, bic = bic, ebic = ebic, sic=sic, do.call(rbind, lapply(x, coefficients)))
    d <- d[order(d[, "k_max"], decreasing = TRUE), , drop = FALSE]
    d <- as.data.frame(d)
    writeLines("HGLM Fit Sequence:")
    print(d, ...)
    invisible(d)
}


connames <- function(x) sapply(x, function(obj) class(obj)[1L])


# constraints <- c(holiglm:::big_m(10), rho_max())
# get_big_m(constraints)
get_big_m <- function(constraints, default = 10) {
    b <- "big_m" == connames(constraints)
    if (any(b)) constraints[[which(b)]][["big_m"]] else default
}


update_big_m <- function(constraints, big_m) {
    b <- "big_m" == connames(constraints)
    if (any(b)) {
        constraints[[which(b)]][["big_m"]] <- big_m
    }
    constraints
}


# constraints = c(rho_max(0.8))
add_constraints <- function(model, constraints, bigM) {
    if (length(constraints) == 0L) return(model)
    bigM_required <- is_bigM_required(constraints)
    if (!any("big_m" == connames(model[["constraints"]])) && any(bigM_required)) {
        if (!any("big_m" == connames(constraints))) {
            constraints <- c(big_m = as.lohglmc(big_m(bigM)), constraints)
        }
        k <- which("big_m" == connames(constraints))
        if (isTRUE(k != 1L)) {
            reorder <- c(k, setdiff(seq_along(constraints), k))
            constraints <- constraints[reorder]
        }
    } else if (!any(bigM_required)) {
        constraints <- constraints[connames(constraints) != "big_m"]
    }

    for (i in seq_along(constraints)) {
        model <- add_constraint(model, constraints[[i]])
    }

    model
}


#' @title Fitting Holistic Generalized Linear Models
#'
#' @description Fit a generalized linear model under constraints.
#'
#' @param model a 'HGLM' model (object of class \code{"hglm_model"}).
#' @param constraints a list of 'HGLM' constraints stored in a list of class \code{"lohglmc"}.
#' @param big_m an upper bound for the coefficients, needed for the big-M constraint.
#'  Required to inherit from \code{"hglmc"}. Currently constraints created by
#'  \code{group_sparsity()}, \code{group_inout()},
#'  \code{include()} and \code{group_equal()} use the big-M set here.
#' @param solver a character string giving the name of the solver to be used for the estimation.
#' @param control a list of control parameters passed to \code{ROI_solve}.
#' @param dry_run a logical; if \code{TRUE} the model is not fit but only constructed.
#' @param approx a logical; if \code{TRUE} uses linear approximation of log-likelihood.
#' @param object_size a character string giving the object size, allowed values
#'  are \code{"normal"} and \code{"big"}. If \code{"big"} is choosen, also
#'  the \pkg{ROI} solution and the \code{"hglm_model"} object are returned.
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' x <- model.matrix(y ~ ., data = dat)
#' model <- hglm_model(x, y = dat[["y"]])
#' fit <- hglm_fit(model, constraints = k_max(3))
#' @return an object of class \code{"hglm.fit"} inheriting from \code{"glm"}.
#' @name hglm_fit
#' @export
hglm_fit <- function(model, constraints = NULL, big_m, solver = "auto", control = list(),
                     dry_run = FALSE, approx = FALSE, object_size = c("normal", "big")) {
    object_size <- match.arg(object_size)
    constraints <- as.lohglmc(constraints)
    if (missing(big_m)) {
        big_m <- get_big_m(constraints, default = 100)
    } else {
        constraints <- update_big_m(constraints, big_m)
    }
    model <- add_constraints(model, constraints, big_m)
    if (isTRUE(dry_run)) return(model)
    tstart <- Sys.time()
    op <- as.OP(x = model)
    op_build_time <- Sys.time() - tstart
    is_mip <- any("big_m" == connames(model[["constraints"]]))
    is_conic <- ROI::is.C_constraint(ROI::constraints(op))
    solver <- choose_solver(solver, family = model$family, is_mip = is_mip, is_conic = is_conic)
    ROI_require_solver(solver)
    # solver = "auto"; control = list()
    tstart <- Sys.time()
    solu <- ROI::ROI_solve(op, solver = solver, control = control)
    op_solve_time <- Sys.time() - tstart
    bicon <- binding_constraints(op, solu)
    if (any(bicon)) {
      L <- L_binding_constraints(op, solu)
      L <- L[, seq_len(NCOL(model$x)), drop = FALSE]
      model$correction_se <- MASS::Null(t(L))
      warning("In hglm_fit: Binding linear constraints detected. ",
              "The standard errors are corrected as described in the vignettes.",
              call. = FALSE)
    }
    fit <- new_hglm_fit(model, roi_solution = solu, object_size = object_size, boundary = any(bicon), approx = approx)
    fit$timing <- c(op_build_time = op_build_time, solve_time = op_solve_time)
    if (any(abs(fit[["coefficients_scaled"]]) >= (big_m - 1e-4)) && is_mip) {
        warning("In hglm_fit: At least one of the Big-M constraints is binding! ",
                "Increase the 'big_m' value and repeat the estimation.")
        fit[["coefficients"]] <- rep.int(NA_real_, length(fit[["coefficients"]]))
    }
    if (!fit$converged) {
        status <- fit[["solution_status"]]
        warning(sprintf("In hglm_fit: Solver '%s' reported the following error: %s",
                        status[["msg"]][["solver"]], status[["msg"]][["message"]]))
    }
    fit
}


# create a new hglm fit object
# currently we use the refit version
# roi_solution <- solu
new_hglm_fit <- function(model, roi_solution, object_size, boundary = FALSE, approx = FALSE) {
    intercept <- has_intercept(model)
    x_scaled <- model$x; x <- unscale_model_data(x_scaled)
    y_scaled <- model$y; y <- unscale_response(y_scaled)
    family <- model$family
    weights <- model$weights; offset <- model$offset; ynames <- model$ynames

    conv <- solution(roi_solution, "status_code") == 0L
    opvals <- solution(roi_solution, force = TRUE)

    idx_bigM <- model[["variable_categories"]][["constr.big_m"]]
    idx_coef <- model[["variable_categories"]][["logli.coef"]]

    if (length(idx_bigM)) {
        if (intercept) {
            coef_indicators <- as.integer(c(1L, opvals[idx_bigM]))
        } else {
            coef_indicators <- as.integer(opvals[idx_bigM])
        }
        coeffi_scaled <- coef_indicators * opvals[idx_coef]  # Use the binary indicators to set to zero.
    } else {
        coef_indicators <- rep.int(1L, NCOL(x))
        coeffi_scaled <- opvals[idx_coef]
    }
    coeffi <- setNames(coeffi_scaled / get_scale(model), colnames(x))

    if (intercept) {
        if (is_probit(family, approx)) {
            correction <- sum(coeffi_scaled[-1] / get_scale(model)[-1] * get_center(model)[-1]) / (2 * sqrt(2 / pi))
        } else {
            correction <- sum(coeffi_scaled[-1] / get_scale(model)[-1] * get_center(model)[-1])
        }
        coeffi[1] <- coeffi_scaled[1] - correction
    }

    # change coefficients due to scaled response
    ym <- get_center_response(model); ys <- get_scale_response(model)
    if (ym != 0 || ys != 1) {
        if (family[["link"]] != "log") {
            coeffi <- coeffi * family$linkfun(ys)
            if (intercept) coeffi[1] <- coeffi[1] + family$linkfun(ym)
        } else {
            if (intercept) coeffi[1] <- coeffi[1] + family$linkfun(ys)
            # else coeffi <- coeffi + abs(sign(round(coeffi, 5))) * family$linkfun(ys) # deactivated scaling response for no intercept
        }
    }

    linkinv <- family$linkinv
    nobs <- NROW(y)
    nvars <- ncol(x)
    is_active <- setNames(coef_indicators != 0, colnames(x))
    n_active_vars <- sum(is_active)
    eta <- setNames(drop(x[, is_active, drop = FALSE] %*% coeffi[is_active]), ynames)
    mu <- setNames(linkinv(eta), ynames)
    mu.eta.val <- family$mu.eta(eta)
    good <- (weights > 0) & (mu.eta.val != 0)
    residuals <- setNames((y - mu) / family$mu.eta(eta), ynames)
    dev <- sum(family$dev.resids(y, mu, weights))

    z <- (eta - offset)[good] + (y - mu)[good] / mu.eta.val[good]
    fam_var <- family$variance(mu)
    # correct zero variance to avoid divison through 0 in w
    fam_var[abs(fam_var) < 1e-08] <- 1e-08
    w <- sqrt((weights[good] * mu.eta.val[good]^2) / fam_var[good])
    wt <- setNames(rep.int(0, nobs), ynames)
    wt[good] <- w^2
    wtdmu <- if (intercept) sum(weights * y) / sum(weights) else linkinv(offset)
    nulldev <- sum(family$dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)

    iter <- niter(roi_solution)
    qr_tolerance <- 1e-07  # glm uses # min(1e-07, control$epsilon / 1000)
    qr <- qr.default(x[, is_active, drop = FALSE] * w, qr_tolerance, LAPACK = FALSE)
    effects <- qr.qty(qr, z * w)
    qr$tol <- qr_tolerance

    nr <- min(sum(good), nvars)
    if (nr < nvars) {
        Rmat <- diag(nvars)
        Rmat[1L:nr, 1L:nvars] <- qr$qr[1L:nr, 1L:nvars]
    } else {
        Rmat <- qr$qr[seq_len(n_active_vars), seq_len(n_active_vars)]
    }
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    dimnames(Rmat) <- list(colnames(qr$qr), colnames(qr$qr))

    R <- Rmat
    rank <- sum(coef_indicators)

    aic.model <- family$aic(y, model$n, mu, weights, dev) + 2 * rank
    if(is.nan(aic.model) || is.na(aic.model)) {
        aic.model <- dev + 2*rank
    }

    resdf <- n.ok - rank

    solution_status <- ROI::solution(roi_solution, "status")

    if (object_size == "big") {
        fit <- list(coefficients = coeffi, coefficients_scaled = coeffi_scaled,
                    residuals = residuals, fitted.values = mu,
                    effects = effects, R = R, rank = rank, qr = qr,
                    family = family, linear.predictors = eta, deviance = dev,
                    aic = aic.model, null.deviance = nulldev, iter = iter,
                    weights = wt, prior.weights = setNames(weights, ynames),
                    df.residual = resdf, df.null = nulldf, y = setNames(y, ynames),
                    converged = conv, boundary = boundary, coefficients.selected = is_active,
                    roi_solution = roi_solution, solution_status = solution_status,
                    hglm_model = model)
    } else {
        fit <- list(coefficients = coeffi, coefficients_scaled = coeffi_scaled,
                    residuals = residuals, fitted.values = mu,
                    effects = effects, R = R, rank = rank, qr = qr,
                    family = family, linear.predictors = eta, deviance = dev,
                    aic = aic.model, null.deviance = nulldev, iter = iter,
                    weights = wt, prior.weights = setNames(weights, ynames),
                    df.residual = resdf, df.null = nulldf, y = setNames(y, ynames),
                    converged = conv, boundary = boundary, coefficients.selected = is_active,
                    solution_status = solution_status)
    }

    structure(fit, class = c("hglm.fit", "glm", "lm"))
}


# variable categories
# - "logli.coef"
# - "logli.aux"
# - "constr.big_m"
# - "constr.equal_sign"

#  returns an object of hglm
#' @title Create a HGLM Model
#'
#' @description Create a HGLM model object.
#'
#' @details No standardization prior to fitting the model takes place. If a x or y standardization is wanted, the user has to do this beforehand.
#'
#' @param x a numeric matrix giving the design matrix.
#' @param y a vector giving the response variables.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param weights an optional vector of 'prior weights' to be used for the estimation.
#' @param frame an optional model frame object.
#' @param approx a logical; if \code{TRUE} uses linear approximation of log-likelihood.
#' @param solver a character string giving the name of the solver to be used for the estimation.
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' x <- model.matrix(y ~ ., data = dat)
#' hglm_model(x, y = dat[["y"]])
#' @return An object of class \code{"hglm_model"}.
#' @name hglm_model
#' @export
hglm_model <- function(x, y, family = gaussian(), weights = NULL, frame = NULL, solver = "auto", approx = FALSE) {
    assert(check_class(family, "family"))
    assert_true(NROW(x) == NROW(y))
    assert_true(NROW(x) > 0)
    approx_name <- if (approx) "_approx" else ""
    # glm.fit like - init - START
    x <- as.matrix(x)
    if (!all(is.finite(x))) {
        col_sum <- colSums(!is.finite(x))
        stop("model.matrix contains non-finite elements in the following columns: \n       ",
             deparse(names(col_sum)[col_sum > 0L]))
    }
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    nobs <- NROW(y)
    nvars <- ncol(x)  # used in in initialize
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    offset <- rep.int(0, nobs)
    variance <- family$variance  # used in in initialize
    linkinv <- family$linkinv  # used in in initialize
    etastart <- NULL  # used in in initialize
    n <- NULL  # fix so R CMD check does not compare.
    eval(family$initialize)
    # glm.fit like - init - END
    assert_numeric(y)
    construct_optimization_problem <- loglike_function(family, approx_name)
    if (is.null(construct_optimization_problem)) {
        msg <- sprintf("hglm model of family '%s' with link '%s' is not implemented!",
                       family$family, family$link)
        stop(msg)
    }
    if (is_power(family)) {
        op <- construct_optimization_problem(x, y, extract_power(family))
    } else {
        op <- construct_optimization_problem(x, y, weights = weights, solver = solver)
    }
    n_aux_variables <- op[["n_of_variables"]] - ncol(x)

    variable_categories <- list("intercept" = which(colnames(x) == "(Intercept)"), 
        "logli.coef" = seq_len(ncol(x)), "logli.aux" = ncol(x) + seq_len(n_aux_variables))
    variable_names <- c(colnames(x), sprintf("%s%i", rep.int("logli.aux", n_aux_variables), seq_len(n_aux_variables)))
    model <- list(x = x, y = y, family = family, weights = weights, frame = frame, n = n,
                  offset = offset, ynames = ynames, nvars = op[["n_of_variables"]],
                  variable_categories = variable_categories, variable_names = variable_names,
                  loglikelihood = op, constraints = list(), types = dense_types(op), scaler = "off")
    class(model) <- "hglm_model"
    model
}


dense_types <- function(x) {
    if (is.null(xty <- types(x))) rep.int("C", x[["n_of_variables"]]) else xty
}


update_nvars <- function(x, nvar) UseMethod("update_nvars")
update_nvars.NULL <- function(x, nvar) x
update_nvars.NO_constraint <- function(x, nvar) x


update_nvars.L_constraint <- function(x, nvar) {
    x[["L"]][["ncol"]] <- nvar
    x[["L"]][["dimnames"]] <- list(NULL, NULL)
    x
}

update_nvars.C_constraint <- function(x, nvar) {
    x[["L"]][["ncol"]] <- nvar
    x[["L"]][["dimnames"]] <- list(NULL, NULL)
    x
}

is_bigM_required <- function(constraints) {
    is_required <- function(x) {
        big_m_needed <- c("k_max", "rho_max", "group_sparsity", "group_inout", "include")
        isTRUE(any(class(x)[1] == big_m_needed)) || (isTRUE(class(x)[1] == "linear") && isTRUE(x[["on_big_m"]]))
    }
    as.logical(lapply(constraints, is_required))
}


#' @title Convert to OP
#'
#' @description Convert an object of class \code{"hglm_model"} into a \pkg{ROI} optimization problem (\code{\link[ROI]{OP}}).
#' @details
#' This function is mainly for internal use and advanced users which want of
#' alter the model object or the underlying optimization problem.
#' This function converts the model object created by \code{\link{hglm_model}}
#' into a conic optimization problem solveable via \code{\link[ROI]{ROI_solve}}.
#' @param x an object inheriting from \code{"hglm_model"}.
#' @return A \pkg{ROI} object of class \code{"OP"}.
#' @examples
#' dat <- rhglm(100, c(1, 2, -3, 4, 5, -6))
#' # Use hglm with option dry_run
#' model <- hglm(y ~ ., data = dat, dry_run = TRUE)
#' op <- as.OP(model)
#' # User hglm_model
#' x <- model.matrix(y ~ ., data = dat)
#' model <- hglm_model(x, dat[["y"]])
#' op <- as.OP(model)
#' @name as.OP.hglm_model
#' @export
as.OP.hglm_model <- function(x) {
    op <- x$loglikelihood
    row_names <- vector("list", length(x[["constraints"]]) + 1L)
    row_names[[1]] <- sprintf("logli_%i", seq_len(op[["n_of_constraints"]]))
    # no additional constraints found, return likelikhood
    if (length(x$constraints) == 0L) {
        op[["row_names"]] <- row_names[[1L]]
        return(op)
    }

    nvars <- max(c(nvar(op), as.integer(lapply(x[["constraints"]], nvar.constraint))))
    for (i in seq_along(x[["constraints"]])) {
        prefix <- names(x[["constraints"]])[i]
        iseq <- seq_len(nrow(x[["constraints"]][[i]]))
        row_names[[i + 1]] <- sprintf("%s_%i", prefix, iseq)
        x[["constraints"]][[i]] <- update_nvars(x[["constraints"]][[i]], nvars)
    }
    constr <- do.call(rbind, x$constraints)
    row_names <- do.call(c, row_names)

    if (!is.null(op$objective$L)) {
        op$objective$L$ncol <- nvars
    }
    # Q needs to be PSD
    if (!is.null(op$objective$Q)) {
        op$objective$Q$nrow <- nvars
        op$objective$Q$ncol <- nvars
    }
    # set length of ROI objective
    attr(op$objective, "nobj") <- nvars
    op$objective$names <- NULL
    op$n_of_variables <- nvars
    op$types <- x$types
    op$constraints$L$ncol <- nvars
    dimnames(op$constraints$L) <- NULL
    dimnames(op$objective$L) <- NULL
    dimnames(op$objective$Q) <- NULL
    dimnames(constr$L) <- NULL
    if (inherits(op$constraints, "NO_constraint")) {
        op$constraints <- constr
    } else {
        op$constraints <- rbind(op$constraints, constr)
    }
    op[["n_of_constraints"]] <- NROW(op$constraints)

    n_new <- nvars - length(x[["variable_names"]])
    variable_names <- c(x[["variable_names"]], sprintf("AUX_%i", seq_len(n_new)))
    op[["objective"]][["names"]] <- variable_names
    op[["constraints"]][["names"]] <- variable_names
    if (op[["n_of_constraints"]] == length(row_names)) {
        op[["row_names"]] <- row_names
    } else {
        stop("dimension missmatch in as.OP")
    }
    op
}


print.hglm_model <- function(x, ...) {
    writeLines(sprintf("Holistic Generalized Linear Model"))
    writeLines(sprintf("  - Family: %s", x$family$family))
    writeLines(sprintf("  - Link function: %s", x$family$link))
    writeLines("  - Constraints")
    if (length(x$constraints) == 0L) {
        writeLines("    + No constaints added.")
    } else {

    }
}


# Creates a new Unique Constraint name
#
# param names a character vector giving the current names.
# param prefix a character string giving the desired prefix.
# 
# return A character string giving the new name.
#
# examples
# constraint_name(names(model[["constraints"]]), "group_sparsity")
constraint_name <- function(names, prefix) {
    ncon <- if (length(names) == 0L) 0L else sum(startsWith(names, prefix))
    paste(prefix, sprintf("%i", ncon + 1L), sep = "_")
}


##
##  p = k + 1
##
# @title Add Big M-Constraint to Model
#
# @description tbd
# @param model a \code{"hglm_model"}, to which the Big M-constraint
#              will be added.
# @param big_m an upper bound for the coefficients, needed for the big-M constraint.
# @return An object of class \code{"hglm_model"}.
# export
add_bigm_constraint <- function(model, big_m = 1000) {
    if (any(names(model[["variable_categories"]]) == "constr.big_m")) {
        return(model)
    }
    p <- ncol(model$x)  # number of covariates
    nobj <- model[["nvars"]]
    vcat <- model[["variable_categories"]]
    io <- as.integer(has_intercept(model))  # intercept offset
    k <- p - io
    j <- seq_len(k)
    big_m <- rep.int(-big_m, k)
    LB <- cbind(stm(j, j + io, rep.int(-1, k)), stzm(k, nobj - p), stdm(big_m))
    UB <- cbind(stm(j, j + io, rep.int( 1, k)), stzm(k, nobj - p), stdm(big_m)) # nolint
    m <- 2 * k
    constr <- L_constraint(rbind(LB, UB), dir = leq(m), rhs = double(m))
    model[["constraints"]][["bigM"]] <- constr
    model[["types"]] <- c(rep("C", nobj), rep("B", k))
    col_idx <- model[["variable_categories"]][["constr.big_m"]]
    model[["variable_categories"]][["constr.big_m"]] <- c(col_idx, nobj + seq_len(k))
    model[["variable_names"]] <- c(model[["variable_names"]], sprintf("BM_%s", setdiff(colnames(model$x), "(Intercept)")))
    model[["nvars"]] <- ncol(constr)
    model
}


# Group Sparsity Constraint
#
# @description Non-exported Helper Function for creatign group sparsity constraints.
# @param model a \code{"hglm_model"}, for which the group sparsity constraint
#              will be defined.
# @param group a vector specifying the group (indices) to which the
#              constraint shall be applied.
# @param k_max an integer giving the maximum number of covariates to be used
#              for the specified group.
# @return An object of class \code{"L_constraint"}.
group_sparsity_constraint <- function(model, group, k_max = 1) {
    stopifnot(all(group %in% model[["variable_categories"]][["logli.coef"]]))
    ones <- rep.int(1, length(group))
    idx <- model[["variable_categories"]][["constr.big_m"]]
    if (is.null(idx)) stop("'bigm_constraint' has to be added before group sparsity constraint!")
    idx <- idx[group - has_intercept(model)]
    L <- stm(ones, idx, ones, ncol = model[["nvars"]])
    L_constraint(L, dir = leq(1L), rhs = k_max + 0.5)
}


# @title Add Group Sparsity Constraint to Model
#
# @description tbd
# @param model a \code{"hglm_model"}, to which the group sparsity constraint
#              will be added.
# @param vars a vector specifying the group (indices) to which the
#              constraint shall be applied.
# @param k_max an integer giving the maximum number of covariates to be used
#              for the specified group.
# @return An object of class \code{"hglm_model"}.
# export
add_group_sparsity_constraint <- function(model, vars, k_max = 1) {
    group <- match_vars(model, vars)[["mm_col_idx"]]
    cname <- constraint_name(names(model[["constraints"]]), "group_sparsity")
    model[["constraints"]][[cname]] <- group_sparsity_constraint(model, group, k_max)
    model
}


# @title Add Global Sparsity Constraint to Model
#
# @description tbd
# @param model a \code{"hglm_model"}, to which the global sparsity constraint
#              will be added.
# @param k_max an integer giving the maximum number of covariates to be used.
# @return An object of class \code{"hglm_model"}.
# export
add_global_sparsity_constraint <- function(model, k_max) {
    group <- seq.int(1 + has_intercept(model), ncol(model$x))
    cname <- constraint_name(names(model[["constraints"]]), "global_sparsity")
    model[["constraints"]][[cname]] <- group_sparsity_constraint(model, group, k_max)
    model
}


# @title Add In-Out Constraint to Model
#
# @description Forces that all covariates in the specified group are
#              either zero or all are nonzero.
# @param model a \code{"hglm_model"}, to which the in-out constraint
#              will be added.
# @param vars a vector specifying the group (indices) to which the
#              constraint shall be applied.
# @return An object of class \code{"hglm_model"}.
# export
add_group_inout_constraint <- function(model, vars) {
    group <- match_vars(model, vars)[["mm_col_idx"]]
    idx <- model[["variable_categories"]][["constr.big_m"]]
    if (is.null(idx)) stop("'bigm_constraint' has to be added before group sparsity constraint!")
    idx <- idx[group - has_intercept(model)]
    i <- rep.int(seq_len(length(group) - 1L), 2)
    j <- c(head(idx, -1), tail(idx, -1))
    v <- rep.int(1, length(group) - 1L)
    L <- stm(i, j, c(v, -v), ncol = model[["nvars"]])
    cname <- constraint_name(names(model[["constraints"]]), "group_inout")
    model[["constraints"]][[cname]] <- L_constraint(L, dir = eq(NROW(L)), rhs = double(NROW(L)))
    model
}


# @title Add In Constraint to Model
#
# @description Ensure that all covariates specified by idx are
#              nonzero.
# @param model a \code{"hglm_model"}, to which the in constraint
#              will be added.
# @param vars an integer vector specifying the indices for covariates which
#            have to be in the model.
# @return An object of class \code{"hglm_model"}.
# export
add_in_constraint <- function(model, vars) {
    indices <- match_vars(model, vars)[["mm_col_idx"]]
    idx <- model[["variable_categories"]][["constr.big_m"]]
    if (is.null(idx)) stop("'bigm_constraint' has to be added before group sparsity constraint!")
    idx <- idx[indices - has_intercept(model)]
    L <- stm(seq_len(length(idx)), idx, rep.int(1, length(idx)), ncol = model[["nvars"]])
    cname <- constraint_name(names(model[["constraints"]]), "idx_in")
    model[["constraints"]][[cname]] <- L_constraint(L, dir = eq(NROW(L)), rhs = rep.int(1, NROW(L)))
    model
}


# @title Add Group Equal Constraint to Model
#
# @description Ensure that all covariates in the specified group have
#              the same value.
# @param model a \code{"hglm_model"}, to which the group_equal constraint
#              will be added.
# @param vars a vector specifying the group (indices) to which the
#              constraint shall be applied.
# @return An object of class \code{"hglm_model"}.
# export
add_group_equal_constraint <- function(model, vars) {
    cbn <- combn(vars, 2L)
    col_names <- sort(unique(as.vector(cbn)))
    j <- match(as.vector(cbn), col_names)
    i <- rep(seq_len(NCOL(cbn)), each = 2L)
    v <- rep.int(c(1, -1), NCOL(cbn))
    L <- as.matrix(slam::simple_triplet_matrix(i, j, v))
    colnames(L) <- col_names
    constr <- .linear_constraint(model, L, dir = rep.int("==", NROW(L)), double(NROW(L)))
    cname <- constraint_name(names(model[["constraints"]]), "group_equal")
    model[["constraints"]][[cname]] <- constr
    model
}


delme__add_group_equal_constraint <- function(model, vars) {
    group <- match_vars(model, vars)[["mm_col_idx"]]
    idx <- model[["variable_categories"]][["constr.big_m"]]
    if (is.null(idx)) stop("'bigm_constraint' has to be added before group sparsity constraint!")
    idx <- idx[group - has_intercept(model)]
    i <- rep.int(seq_len(length(group) - 1L), 2)
    j <- c(head(idx, -1), tail(idx, -1))
    v <- rep.int(1, length(group) - 1L)
    L <- stm(i, j, c(v, -v), ncol = model[["nvars"]])
    constr <- L_constraint(L, dir = eq(NROW(L)), rhs = double(NROW(L)))
    cname <- constraint_name(names(model[["constraints"]]), "group_equal")
    model[["constraints"]][[cname]] <- constr
    model
}


add_bound <- function(model, kvars, btype = c("lower", "upper")) {
    btype <- match.arg(btype)
    mapping <- match_kvars(model, kvars)
    vals <- mapping[["value"]]
    idx <- mapping[["mm_col_idx"]]
    L <- simple_triplet_matrix(i = seq_along(idx), j = idx, v = rep.int(1, length(idx)),
                               ncol = model[["nvars"]])
    xeq <- if (btype == "lower") ROI::geq else ROI::leq
    constr <- L_constraint(L, dir = xeq(length(idx)), rhs = get_scale(model)[idx] * vals / get_scale_response(model))
    bound_type <- sprintf("%s_bound", btype)
    cname <- constraint_name(names(model[["constraints"]]), bound_type)
    model[["constraints"]][[cname]] <- constr
    model
}


# @title Add Lower Bound
#
# @description 
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param kvars TODO
# @return An object of class \code{"hglm_model"}.
# export
add_lower_bound <- function(model, kvars) {
    add_bound(model, kvars, "lower")
}


# @title Add Upper Bound
#
# @description 
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param kvars TODO
# @return An object of class \code{"hglm_model"}.
# export
add_upper_bound <- function(model, kvars) {
    add_bound(model, kvars, "upper")
}


# This is not perfect since it is easier for the user this way but, the
# user can not combine linear constraints on the integer and other variables.
.linear_constraint <- function(model, kvars, dir, rhs, on_big_m = FALSE) {
    vars <- colnames(kvars)
    mapping <- match_vars(model, vars)
    idx <- mapping[["mm_col_idx"]]
    name_mapping <- mapping[["vars"]]
    if (on_big_m) {
        # for the big M constraint we need no scaling
        mat <- if (is.null(name_mapping)) kvars else kvars[, name_mapping, drop = FALSE]
    } else {
        diag_values <- get_scale_response(model) / get_scale(model)[idx]
        diag_mat <- diag(diag_values, nrow = length(diag_values), ncol = length(diag_values))
        if (is.null(name_mapping)) {
            mat <- kvars %*% diag_mat
        } else {
            mat <- kvars[, name_mapping, drop = FALSE] %*% diag_mat
        }
    }
    p <- ncol(model$x)  # number of covariates
    k <- p - has_intercept(model)
    nobj <- ncol(constraints(model$loglikelihood))
    L <- as.simple_triplet_matrix(mat)
    L[["ncol"]] <- nobj + k
    if (on_big_m) {
        j <- idx[L[["j"]]]
        if (length(model[["variable_categories"]][["intercept"]]) > 0) {
            L[["j"]] <- model[["variable_categories"]][["constr.big_m"]][j - 1L]
        } else {
            L[["j"]] <- model[["variable_categories"]][["constr.big_m"]][j]
        }
    } else {
        L[["j"]] <- idx[L[["j"]]]
    }
    return(L_constraint(L, dir = dir, rhs = rhs))
}


# @title Add Linear Constraint
#
# @description
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param kvars TODO
# @param dir TODO
# @param rhs TODO
# @return An object of class \code{"hglm_model"}.
# @examples
# \dontrun{
# add_linear_constraint(model, c("a" = 3, "b" = 4, "c" = 6), ">=", 3)
# }
# export
add_linear_constraint <- function(model, L, dir, rhs, on_big_m = FALSE) {
    if (is.vector(L)) {
        L <- matrix(L, nrow = 1L, dimnames = list(NULL, names(L)))
    }
    if (length(dir) == 1L && NROW(L) > 1L) {
        dir <- rep.int(dir, NROW(L))
        rhs <- rep.int(rhs, NROW(L))
    }
    assert_true(NROW(L) == length(rhs))
    if (on_big_m) {
        cname <- constraint_name(names(model[["constraints"]]), "linear_constraint_bigM")
    } else {
        cname <- constraint_name(names(model[["constraints"]]), "linear_constraint")
    }
    model[["constraints"]][[cname]] <- .linear_constraint(model, kvars = L, dir, rhs, on_big_m = on_big_m)
    model
}


# @title Add Equal Sign Constraint
#
# @description
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param kvars TODO
# @param big_m TODO
# @param eps a double giving the epsilon for the equal sign constraint.
#  Since most numerical solver can only handle constraints up to some epsilon,
#  e.g., the constraint \deqn{A x \geq b} typically is only solved to
#  \deqn{|A x - b| \geq 0}.
# @return An object of class \code{"hglm_model"}.
# export
add_sign_constraint <- function(model, kvars, big_m = 100, eps = 1e-6) {
    if (is.character(kvars)) {
        kvars <- setNames(rep.int(1L, length(kvars)), kvars)
    }
    assert_kvars(kvars)
    mapping <- match_kvars(model, kvars)
    idx <- mapping[["mm_col_idx"]]
    nidx <- length(idx)
    i <- seq_along(idx)
    j <- c(idx, rep.int(model[["nvars"]] + 1L, nidx))
    v <- c(mapping[["value"]], rep.int(-big_m, nidx))
    L <- simple_triplet_matrix(c(i, i), j, v = v, ncol = model[["nvars"]] + 1L)
    constr <- L_constraint(L = rbind(L, L), dir = c(leq(nidx), geq(nidx)),
                           rhs = c(rep.int(-eps, nidx), rep.int(-big_m + eps, nidx)))
    cname <- constraint_name(names(model[["constraints"]]), "equal_sign")
    model[["constraints"]][[cname]] <- constr
    model[["types"]] <- c(model[["types"]], "B")
    col_idx <- model[["variable_categories"]][["constr.equal_sign"]]
    model[["variable_categories"]][["constr.equal_sign"]] <- c(col_idx, model[["nvars"]] + 1L)
    vnam <- model[["variable_names"]]
    model[["variable_names"]] <- c(vnam, sprintf("SC_%s", paste(vnam[idx], collapse = ":")))
    model[["nvars"]] <- model[["nvars"]] + 1L
    model
}


# @title Add Pairwise Multicollinearity Constraint to Model
#
# @description Ensure that model is free from extreme pairwise
#              multicollinearity.
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param rho_max a value in the range [0,1] specifying, the maximum
#              allowed collinearity between pairs of covariates
# @param exclude TODO
# @param use an optional character string giving a method for computing
#            covariances in the presence of missing values.
#            The parameter is passed to \code{\link[stats]{cor}},
#            therefore see \code{\link[stats]{cor}} for more information.
# @param method a character string indicating which correlation coefficient
#            is to be computed. See \code{\link[stats]{cor}} for more information.
#            The parameter is passed to \code{\link[stats]{cor}},
#            therefore see \code{\link[stats]{cor}} for more information.
# @return An object of class \code{"hglm_model"}.
# export
add_pairwise_multicollinearity_constraint <- function(
    model,
    rho_max = 0.8,
    exclude = "(Intercept)",
    use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"),
    method = c("pearson", "kendall", "spearman")) {
    use <- match.arg(use)
    method <- match.arg(method)
    idx <- match(exclude, colnames(model[["x"]]))
    stopifnot(!anyNA(idx))
    rho <- cor(model$x[ , -idx, drop = FALSE])
    M <- abs(rho ) > rho_max
    df <- as.data.frame(unclass(as.simple_triplet_matrix(M))[1:3])
    df <- df[df[["i"]] > df[["j"]], ]
    if (NROW(df)) {
        mapping_rho_coef <- match(colnames(rho), setdiff(colnames(model$x), "(Intercept)"))
        mapping_rho_zvar <- model[["variable_categories"]][["constr.big_m"]][mapping_rho_coef]
        i <- rep.int(seq_len(NROW(df)), 2)
        j <- mapping_rho_zvar[c(df[["i"]], df[["j"]])]
        v <- rep.int(1, 2 * NROW(df))
        L <- simple_triplet_matrix(i, j, v, ncol = model[["nvars"]])
        constr <- L_constraint(L, dir = leq(NROW(L)), rhs = rep.int(1, NROW(L)))
    } else {
        constr <- NULL
    }
    cname <- constraint_name(names(model[["constraints"]]), "pairwise_multicollinearity")
    model[["constraints"]][[cname]] <- constr
    model
}


# @title Add Pairwise Sign Constraint to Model
#
# @description Ensure that variables which exhibit strong pairwise correlation
#  have the same sign.
#
# @param model a \code{"hglm_model"}, to which the pairwise
#              multicollinearity constraint will be added.
# @param rho_max a value in the range [0,1] specifying, the maximum
#              allowed collinearity between pairs of covariates
# @param exclude a character vector giving the names of variables to be excluded
#        from the correlation calculations.
# @param big_m TODO
# @param eps a double giving epsilon to be added 
# @param use an optional character string giving a method for computing
#            covariances in the presence of missing values.
#            The parameter is passed to \code{\link[stats]{cor}},
#            therefore see \code{\link[stats]{cor}} for more information.
# @param method a character string indicating which correlation coefficient
#            is to be computed. See \code{\link[stats]{cor}} for more information.
#            The parameter is passed to \code{\link[stats]{cor}},
#            therefore see \code{\link[stats]{cor}} for more information.
# @return An object of class \code{"hglm_model"}.
# export
add_pairwise_multicollinearity_equal_sign_constraint <- function(
    model,
    rho_max = 0.8,
    exclude = "(Intercept)",
    big_m = 100,
    eps = 1e-6,
    use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"),
    method = c("pearson", "kendall", "spearman")) {
    use <- match.arg(use)
    method <- match.arg(method)
    idx <- match(exclude, colnames(model[["x"]]))
    stopifnot(!anyNA(idx))
    X <- model$x[ , -idx, drop = FALSE]
    rho <- cor(X, use = use, method = method)
    M <- abs(rho) > rho_max
    df <- as.data.frame(unclass(as.simple_triplet_matrix(M))[1:3])
    df <- df[df$i > df$j,]
    if (NROW(df)) {
        for (k in seq_len(NROW(df))) {
            i <- df[k, "i"]
            j <- df[k, "j"]
            kvars <- setNames(c(1, sign(rho[i, j])), c(colnames(rho)[i], colnames(rho)[j]))
            model <- add_sign_constraint(model, kvars, big_m = big_m, eps = eps)
        }
    }
    model
}


#' Update the Model Object
#'
#' @param model an object inheriting from \code{"hglm_model"}.
#' @param op an \code{\link[ROI]{OP}} giving the new objective.
#'
#' @export
update_objective <- function(model, op) {
    checkmate::assert_class(model, "hglm_model")
    checkmate::assert_class(op, "OP")
    p <- ncol(model$x)
    model[["loglikelihood"]] <- op
    nvars <- model[["nvars"]] <- op[["n_of_variables"]]
    ind_logli_aux <- seq(p + 1L, model[["nvars"]])
    model[["variable_categories"]][["logli.aux"]] <- ind_logli_aux
    model[["variable_names"]] <- c(model[["variable_names"]][seq_len(p)],
                                   sprintf("logli.aux%i", seq_len(nvars - p)))
    model[["types"]] <- dense_types(op)
    model
}

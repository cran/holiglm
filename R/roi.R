
nvar <- function(x) UseMethod("nvar")
nvar <- function(x) x[["n_of_variables"]]

ncon <- function(x) UseMethod("ncon")
ncon <- function(x) x[["n_of_constraints"]]


dim.OP <- function(x) c(ncon(x), nvar(x))


niter <- function(roi_solution) {
    .niter <- function(roi_solution) {
        solver <- solution(roi_solution, "status")[["msg"]][["solver"]]
        msg <- solution(roi_solution, "msg")
        n_iter <- switch(solver,
            "gurobi" = msg[["itercount"]],
            "ecos" = msg[["retcodes"]][["iter"]],
            "mosek" = msg[["iinfo"]][["INTPNT_ITER"]],
            "osqp " = msg[["info"]][["iter"]],
            "quadprog" = sum(msg[["iterations"]]),
            NA_integer_
        )
        return(n_iter)
    }
    tryCatch(.niter(roi_solution), error = function(e) NA_integer_)
}


binding_constraints <- function(op, solu, tol = 1e-4) {
    k <- grep("^(linear_constraint|lower_bound|upper_bound)_\\d+", op[["row_names"]])
    if (length(k) == 0L) return(NULL)
    if (is.C_constraint(constraints(op))) {
        con <- constraints(op)
        L <- con[["L"]][k,]
        dir <- c("==", "<=")[con[["cones"]][["cone"]][k]]
        rhs <- con[["rhs"]][k]
    } else if (is.L_constraint(constraints(op))) {
        con <- constraints(op)
        L <- con[["L"]][k,]
        dir <- con[["dir"]][k]
        rhs <- con[["rhs"]][k]
    } else {
        stop("unsupported constraint type")
    }
    s <- solution(solu, force = TRUE)
    abs(drop(as.matrix(L) %*% s) - rhs) < tol
}


L_binding_constraints <- function(op, solu, tol = 1e-4) {
    k <- grep("^(linear_constraint|group_equal|lower_bound|upper_bound)_\\d+", op[["row_names"]])
    if (length(k) == 0L) return(NULL)
    if (is.C_constraint(constraints(op))) {
        con <- constraints(op)
        L <- as.matrix(con[["L"]][k, ])
        dir <- c("==", "<=")[con[["cones"]][["cone"]][k]]
        rhs <- con[["rhs"]][k]
    } else if (is.L_constraint(constraints(op))) {
        con <- constraints(op)
        L <- as.matrix(con[["L"]][k, ])
        dir <- con[["dir"]][k]
        rhs <- con[["rhs"]][k]
        if (any(b <- ">=" == dir)) {
            L[b, ] <- -L[b, ]
            rhs[b] <- -rhs[b]
        }
    } else {
        stop("unsupported constraint type")
    }
    s <- solution(solu, force = TRUE)
    is_bind <- abs(drop(L %*% s) - rhs) < tol
    L[dir == "<=", ] <- -L[dir == "<=", ]
    L[is_bind, , drop = FALSE]
}

             
ROI_is_registered <- function(solver) {
    isTRUE(solver %in% names(ROI::ROI_registered_solvers()))
}


ROI_is_installed <- function(solver) {
    isTRUE(solver %in% names(ROI::ROI_installed_solvers()))
}


# Require a Solver
#
# The solver will be loaded
#
# solver <- "ecos"
ROI_require_solver <- function(solver) {
    checkmate::assert_character(solver, len = 1L, any.missing = FALSE)
    if (ROI_is_registered(solver)) return(invisible(0))
    # The gsub is needed for e.g., "nloptr.lbfgs"
    plugin_name <- sprintf("ROI.plugin.%s", gsub("\\..*", "", solver))
    if (ROI_is_installed(solver)) {
        requireNamespace(plugin_name, quietly = TRUE)
    } else {
        stop(sprintf("'%s' can not be found among the installed solvers ", plugin_name),
                     "(in `ROI_installed_solvers()`) please make sure that is installed.")
    }
    return(invisible(0))
}


#' @title Extract Solution
#'
#' @description The solution of the underlying optimization problem, 
#'  can be accessed via the method \code{'solution'}.
#'
#' @param x an object of type \code{'hglm'}.
#' @param type a character giving the name of the solution to be extracted.
#' @param force a logical to control the return value in the case that the
#'        status code is equal to 1 (i.e. something went wrong).
#'        By default force is \code{FALSE} and a solution is only provided
#'        if the status code is equal to 0 (i.e. success). If force is \code{TRUE}
#'        the status code is ignored and solutions are returned also
#'        where the solver signaled an issue.
#' @param ... further arguments passed to or from other methods.
#' @return the extracted solution.
#' @export
solution.hglm <- function(x, type = c("primal", "dual", "aux", "psd", "msg", "objval", "status", "status_code"), force = FALSE, ...) {
    ROI::solution(x[["roi_solution"]])
}


#' @export
solution.hglm.fit <- solution.hglm

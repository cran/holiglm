#' @title Holistic Generalized Linear Models Package
#'
#' @description
#' The holistic generalized linear models package simplifies
#' estimating generalized linear models under constraints.
#' The constraints can be used to,
#' \itemize{
#'   \item bound the domains of specific covariates,
#'   \item impose linear constraints on the covariates,
#'   \item induce sparsity via best subset selection,
#'   \item impose sparsity on groups of variables,
#'   \item restrict the pairwise correlation between the selected coefficients,
#'   \item impose sign coherence constraints on selected covariates and
#'   \item force all predictors within a group either to be selected or not.
#' }
#'
#' This sophisticated constraints are internally implemented via conic optimization.
#' However, the package is designed such that the user, is not required to
#' be familiar with conic optimization but is only required to have basic \R knowledge.
#'
#' @author
#' \itemize{
#'   \item Benjamin Schwendinger (\strong{Maintainer} \email{benjaminschwe@gmail.com})
#'   \item Florian Schwendinger
#'   \item Laura Vana
#' }
#'
#' @references
#' \strong{Holistic regression} \cr
#' Schwendinger, B., Schwendinger, F., & Vana, L. (2022).
#' Holistic Generalized Linear Models.
#' \doi{10.48550/ARXIV.2205.15447}.
#'
#' Bertsimas, D., & King, A. (2016).
#' OR Forum-An Algorithmic Approach to Linear Regression
#' Operations Research 64(1):2-16.
#' \doi{10.1287/opre.2015.1436}
#'
#' Bertsimas, D., & Li, M. L. (2020).
#' Scalable Holistic Linear Regression.
#' Operations Research Letters 48 (3): 203–8.
#' \doi{10.1016/j.orl.2020.02.008}.
#'
#' \strong{Constrained regression} \cr
#' McDonald, J. W., & Diamond, I. D. (1990).
#' On the Fitting of Generalized Linear Models with Nonnegativity Parameter Constraints.
#' Biometrics, 46 (1): 201–206.
#' \doi{10.2307/2531643}
#'
#' Slawski, M., & Hein, M. (2013).
#' Non-negative least squares for high-dimensional linear models: Consistency and sparse recovery without regularization.
#' Electronic Journal of Statistics, 7: 3004-3056.
#' \doi{10.1214/13-EJS868}
#'
#' Carrizosa, E., Olivares-Nadal, A. V., & Ramírez-Cobo, P. (2020).
#' Integer Constraints for Enhancing Interpretability in Linear Regression.
#' SORT. Statistics and Operations Research Transactions, 44: 67-98.
#' \doi{10.2436/20.8080.02.95}.
#'
#' Lawson, C. L., & Hanson, R. J. (1995).
#' Solving least squares problems. Society for Industrial and Applied Mathematics.
#' Society for Industrial and Applied Mathematics.
#' \doi{10.1137/1.9781611971217}
#'
#' \strong{Generalized Linear Models} \cr
#' McCullagh, P., & Nelder, J. A. (2019).
#' Generalized Linear Models (2nd ed.)
#' Routledge.
#' \doi{10.1201/9780203753736}.
#'
#' \strong{Conic Optimization} \cr
#' Boyd, S., & Vandenberghe, L. (2004).
#' Convex Optimization (1st ed.)
#' Cambridge University Press.
#' \url{https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf}.
#' \doi{10.1017/cbo9780511804441}
#'
#' Theußl, S., Schwendinger, F., & Hornik, K. (2020).
#' ROI: An Extensible R Optimization Infrastructure.
#' Journal of Statistical Software 94 (15): 1–64.
#' \doi{10.18637/jss.v094.i15}.
#'
#' @seealso
#' \code{\link{hglm}}, \code{\link{holiglm}}
#'
#' @docType package
#' @name holiglm-package
#'
#' @import slam
#' @import stats
#' @import ROI
#' @import ROI.plugin.ecos
#' @import checkmate
#' @importFrom utils head modifyList tail str combn capture.output
#' @aliases holiglm-package
"_PACKAGE"


#' @export
ROI::as.OP


#' @export
ROI::solution

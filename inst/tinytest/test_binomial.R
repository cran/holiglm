# q("no"); R
# library(tinytest)
suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("ROI"))
applicable_solvers <- c("ecos")
installed_solvers <- ROI_installed_solvers()
solver <- sample(installed_solvers[names(installed_solvers) %in% applicable_solvers], 1L)
require(solver, character.only = TRUE)

test_equalish <- function(x, y, eps = 1e-9) max(abs(x - y)) < eps
except_equalish <- function(x, y, eps = 1e-9) expect_true(test_equalish(x, y, eps))

test <- expression({
eps <- 1e-3  # switch(family$link, probit = e-2, 1e-3)

#' ## Test likelihood function
data <- PlantGrowth[as.integer(PlantGrowth$group) > 1L,]
x <- model.matrix(group ~ ., data = data)
y <- as.integer(data$group) - 2L
loglike_fun <- holiglm:::loglike_function(family)
op <- loglike_fun(x, y)
sol <- ROI_solve(op)
coef0 <- unname(head(sol$solution, ncol(x)))
coef0 <- holiglm:::approx_inverse(coef0, family)
start <- if (family$link == "log") c(-1, double(length(coef0) - 1L))  else NULL
coef1 <- coef(glm(group ~ ., data = data, family = family, start = start))
except_equalish(coef0, coef1, eps)
# cbind(coef0, coef1)

set.seed(0)
data <- data.frame(y = rbinom(500, 1, 0.7),
                   a = factor(sample(1:4, 500, TRUE)),
                   b = factor(sample(1:3, 500, TRUE)))
mm <- model.matrix(y ~ ., data = data)
beta <- c(1, -2, 3, -4, 5, -6)
mu <- exp(mm %*% beta) / (exp(mm %*% beta) + 1)

mm <- model.matrix(y ~ ., data = data)
start <- c(-1, double(ncol(mm) - 1))
coef_glm <- coef(glm(y ~ ., data = data, family = family))

#' ### Without weights
fit <- hglm(y ~ ., data = data, family = family, constraints = NULL, scaler = scaler)
expect_equal(coef(fit), coef_glm, tolerance=1e-4)
# max(abs(coef(fit) - coef_glm))
# cbind(coef(fit), coef_glm)

#' ### Test y matrix
dat2 <- agg_binomial(y ~ ., data, as_list = FALSE)
fit <- hglm(cbind(success, failure) ~ ., data = dat2, family = family, constraints = NULL, scaler = scaler, scale_response = FALSE)
expect_equivalent(coef(fit), coef_glm, tolerance=1e-4)

#' ### Test y probabilities with weights
weights <- dat2$failure + dat2$success
fit <- hglm(success / (failure + success) ~ ., data = dat2, family = family, weights = weights, constraints = NULL,
    scaler = scaler, scale_response = FALSE)
expect_equivalent(coef(fit), coef_glm, tolerance=1e-4)
})



## todo: add expect_warning where needed
scaling <- c("off", "center_standardization", "center_minmax", "standardization", "minmax")
families <- list(binomial(link = "log"), binomial(link = "logit"))

for (scaler in scaling) {
    for (family in families) {
        eval(test)
    }
}

#' ### predict
data <- PlantGrowth[as.integer(PlantGrowth$group) > 1L,]
model <- hglm(group ~ ., data = data, family = binomial())
gmodel <- glm(group ~ ., data = data, family = binomial())
expect_equal(predict(model, type="link"), predict(gmodel, type="link"), tolerance=1e-4)
expect_equal(predict(model, newdata=data[6:15,], type="link"), predict(gmodel, newdata=data[6:15,], type="link"), tolerance=1e-4)
expect_equal(predict(model, type="response"), predict(gmodel, type="response"), tolerance=1e-4)
expect_equal(predict(model, newdata=data[6:15,], type="response"), predict(gmodel, newdata=data[6:15,], type="response"), tolerance=1e-4)

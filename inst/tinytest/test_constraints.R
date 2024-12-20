if (interactive()) {
    library("tinytest")
}


suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("ROI"))

applicable_solvers <- c("ecos")
installed_solvers <- ROI_installed_solvers()
solver <- sample(installed_solvers[names(installed_solvers) %in% applicable_solvers], 1L)
require(solver, character.only = TRUE)


# Constraints
# 1. add_group_sparsity_constraint
# 2. add_global_sparsity_constraint
# 3. add_group_inout_constraint
# 4. add_in_constraint
# 5. add_group_equal_constraint
# 6. add_pairwise_multicollinearity_constraint
# 7. add_sign_constraint
# 8. add_pairwise_multicollinearity_equal_sign_constraint

## 1. Constraints - group sparsity constraint
set.seed(0)
beta <- c(1, 2, 0, 4, 5, 0)
dat <- rhglm(100, beta)
# head(dat)
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
fit <- hglm_fit(model, group_sparsity(c("x1", "x2", "x5"), 1L), big_m = 100)
expect_equivalent(coef(fit)[c("x2", "x5")], c(0, 0))


## 2. Constraints - global sparsity constraint
## add global sparsity constraint
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
fit <- hglm_fit(model, k_max(sum(beta != 0) - 1L), big_m = 100)
expect_equivalent(coef(fit) == 0, beta == 0)


# 3. Constraints - group inout constraint
# add_group_inout_constraint
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
constraints <- c(k_max(sum(beta != 0) - 1L), group_inout(c("x1", "x2", "x3")))
fit <- hglm_fit(model, constraints)
expect_equivalent(coef(fit)[c("x1", "x2", "x3")], double(3))


# 4. Constrain - in-constraint
# add_in_constraint
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
constraints <- c(k_max(sum(beta != 0) - 1L), include(c("x1", "x3")))
fit <- hglm_fit(model, constraints)
expect_equivalent(coef(fit)["x2"], 0)
expect_true(all(coef(fit)[c("x1", "x3")] != 0))


# 5. add_group_equal_constraint
# add_group_equal_constraint
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
fit <- hglm_fit(model, group_equal(vars = c("x1", "x3")))
expect_equivalent(coef(fit)["x1"], coef(fit)["x3"], tolerance = 1e-6)


## 6. Constraints - pairwise multicollinearity constraint
## add pairwise multicollinearity constraint
set.seed(0)
beta <- 1:3
Sigma <- cov_matrix(k = length(beta) - 1L, 1, 2, 0.9)
dat <- rhglm(100, beta, sigma = Sigma)
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
fit <- hglm_fit(model, rho_max(0.8), big_m = 100)
is_none_zero <- coef(fit) > 0
expect_equivalent(sum(is_none_zero[-1L]), 1)
expect_equivalent(coef(fit)[is_none_zero], coef(glm(y ~ x2, data = dat)),
                            tolerance = 1e-4)


## 7. Constraints - equal sign constraint
## add_sign_constraint

# same sign
set.seed(0)
beta <- c(1, 0.5, -1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1", "x2"), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, -1, 1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1", "x2"), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, 1, 1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1", "x2"), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, -1, -1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1", "x2"), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(sign(coef(fit))["x1"], sign(coef(fit))["x2"])


# opposite sign
set.seed(0)
beta <- c(1, 0.5, -1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1" = -1, "x2" = 1), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(-sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, -1, 1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1" = 1, "x2" = -1), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(-sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, 1, 1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1" = -1, "x2" = 1), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(-sign(coef(fit))["x1"], sign(coef(fit))["x2"])


set.seed(0)
beta <- c(1, -1, -1)
dat <- rhglm(100, beta)
model <- hglm(y ~ ., data = dat, constraints = NULL, dry_run = TRUE)
fit <- hglm_fit(model, sign_coherence(c("x1" = -1, "x2" = 1), big_m = 10))
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(-sign(coef(fit))["x1"], sign(coef(fit))["x2"])


## 8. 
## add_pairwise_multicollinearity_equal_sign_constraint
set.seed(1)
beta <- c(1, -2, 3, 4, -5)
Sigma <- cov_matrix(k = length(beta) - 1L, c(1, 3), c(2, 4), c(0.9, 0.95))
Sigma[Sigma == 0] <- round(runif(sum(Sigma == 0), 0, 0.3), 2)
Sigma <- (Sigma + t(Sigma)) / 2
dat <- rhglm(100, beta, sigma = Sigma)
model <- hglm(y ~ ., constraints = NULL, data = dat, dry_run = TRUE)
fit <- hglm_fit(model, pairwise_sign_coherence())
cbind(coef(fit), round(coef(fit), 4), sign(coef(fit)))
expect_equivalent(sign(coef(fit))["x1"], sign(coef(fit))["x2"])
expect_equivalent(sign(coef(fit))["x3"], sign(coef(fit))["x4"])


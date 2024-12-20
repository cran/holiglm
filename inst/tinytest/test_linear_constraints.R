if (FALSE) {
    q("no")
    R
    library("tinytest")
}


suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("ROI.plugin.ecos"))


#
# lower bound
#
set.seed(0)
beta <- c(1, -2, 3)
dat <- rhglm(100, beta)
# coef(glm(y ~ ., data = dat))


fit <- hglm(y ~ ., data = dat, constraints = c(linear(c(x1 = 2, x2 = 1), "==", 0), rho_max(1)))
coefs <- coef(fit)
expect_true(abs(2 * coefs["x1"] + coefs["x2"]) < 1e-8)


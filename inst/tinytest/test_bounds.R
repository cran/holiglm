suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
library("checkmate")
suppressPackageStartupMessages(library("ROI.plugin.ecos"))


#
# lower bound
#
set.seed(0)
beta <- c(1, -2, 3, -1, 4)
dat <- rhglm(100, beta)

scaling <- c("off", "center_standardization", "center_minmax", "standardization", "minmax")

for (scaler in scaling) {
    fit <- hglm(y ~ ., data = dat, constraints = lower(c(x1 = 0, x2 = 0, x3 = 0)), scaler = scaler)
    expect_true(all(coef(fit)[c("x1", "x2", "x3")] >= 0))
}


#
# upper bound
#
for (scaler in scaling) {
    fit <- hglm(y ~ ., data = dat, constraints = c(upper(c(x1 = -3, x2 = 2, x4 = 2)), rho_max(1)), scaler=scaler)
    expect_true(all(coef(fit)[c("x1", "x2", "x3")] < (c(-3, 2, 2) + 1e-6 )))
}

if (interactive()) {
    library("tinytest")
}

Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("holiglm"))

set.seed(0)
data <- data.frame(y = rbinom(50, 1, 0.7),
                   a = round(runif(50), 2),
                   b = factor(sample(1:4, 50, TRUE)),
                   c = factor(sample(1:2, 50, TRUE)),
                   d = round(runif(50), 2))


model <- hglm(y ~ . + poly(a, 3), data = data, constraints = NULL, dry_run = TRUE)

mdf <- holiglm:::match_vars(model, c("a", "b"))
expect_equal(mdf[["mm_col_idx"]], 2:5)

expect_error(holiglm:::match_vars(model, c("c", "unkown_column")))


kvars <- c(a = 3, b = 4)
mdf <- holiglm:::match_kvars(model, kvars)
expect_equal(mdf[["mm_col_idx"]], 2:5)
expect_equal(setNames(mdf[["value"]], mdf[["mm_name"]]), c(a = 3, b2 = 4, b3 = 4, b4 = 4))


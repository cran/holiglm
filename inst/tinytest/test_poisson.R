suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("ROI"))
applicable_solvers <- c("ecos")
installed_solvers <- ROI_installed_solvers()
solver <- sample(installed_solvers[names(installed_solvers) %in% applicable_solvers], 1L)
suppressPackageStartupMessages(require(solver, character.only = TRUE))
## Likelihood poisson

## setup data
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
dobson <- data.frame(outcome, treatment, counts)

## poisson log link function
## example 1
family <- poisson(link="log")
x <- model.matrix(gear ~ mpg, data = mtcars)
y <- mtcars$gear
model <- holiglm:::loglike_poisson_log(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(0.98935, 0.01549)
expect_equal(coef0, coef1, tolerance=1e-4)

## example 2
x <- model.matrix(counts ~ ., data = dobson)
y <- dobson$counts
model <- holiglm:::loglike_poisson_log(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(3.04452, -0.45426, -0.29299, 0, 0)
expect_equal(coef0, coef1, tolerance=1e-4)


## poisson identity link function
## example 1
family <- poisson(link="identity")
x <- model.matrix(gear ~ mpg, data = mtcars)
y <- mtcars$gear
model <- holiglm:::loglike_poisson_identity(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(2.47682, 0.06026)
expect_equal(coef0, coef1, tolerance=1e-4)

## example 2
x <- model.matrix(counts ~ ., data = dobson)
y <- dobson$counts
model <- holiglm:::loglike_poisson_identity(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(21.5307, -7.76271, -5.38842, -0.59051, -0.85045)
expect_equal(coef0, coef1, tolerance=1e-4)


## poisson sqrt link function
## example 1
family <- poisson(link="sqrt")
x <- model.matrix(gear ~ mpg, data = mtcars)
y <- mtcars$gear
model <- holiglm:::loglike_poisson_sqrt(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(1.61049, 0.01531)
expect_equal(coef0, coef1, tolerance=1e-4)

## example 2
x <- model.matrix(counts ~ ., data = dobson)
y <- dobson$counts
model <- holiglm:::loglike_poisson_sqrt(x, y)
sol <- ROI_solve(model)
coef0 <- unname(head(sol$solution, ncol(x)))
coef1 <- c(4.61421, -0.93424, -0.62635, -0.03605, -0.05436)
expect_equal(coef0, coef1, tolerance=1e-4)


# ## poisson power link function
# ## example 1
# x <- model.matrix(gear ~ mpg, data = mtcars)
# y <- mtcars$gear
# mus <- c(1/2, 1/3, 1/4, 1/5)
# res_coefs <- structure(c(1.61049, 1.37992, 1.27513, 1.2155, 0.01531, 0.00814, 0.00545, 0.00407), .Dim = c(4L, 2L))
# poisson_power_modell <- function(mu, x, y) {
#     family <- poisson(link=power(mu))
#     model <- holiglm:::loglike_poisson_power(x, y, mu)
#     sol <- ROI_solve(model)
#     unname(head(sol$solution, ncol(x)))
# }

# for (i in seq_along(mus)) {
#     expect_equal(poisson_power_modell(mus[i], x, y), res_coefs[i,], tolerance=1e-4)
# }

# ## example 2
# x <- model.matrix(counts ~ ., data = dobson)
# y <- dobson$counts
# mus <- c(1/2, 1/3, 1/4, 1/5)
# res_coefs <- structure(c(4.61421, 2.76759, 2.14452, 1.84053, -0.93424, -0.38825, -0.23001, -0.15975, -0.62635, -0.25706,
#     -0.15132, -0.10469, -0.03605, -0.00996, -0.00441, -0.00245, -0.05436, -0.01523, -0.00679, -0.00378), .Dim = 4:5)

# for (i in seq_along(mus)) {
#     expect_equal(poisson_power_modell(mus[i], x, y), res_coefs[i,], tolerance=1e-4)
# }

# hglm_model
x <- model.matrix(counts ~ ., data = dobson)
y <- dobson$counts
model <- hglm_model(x, y, family=poisson())
op <- as.OP(model)
sol <- ROI_solve(op)
coef0 <- unname(head(solution(sol), NCOL(x)))
coef1 <- c(3.04452, -0.45426, -0.29299, 0, 0)
expect_equal(coef0, coef1, tolerance=1e-4)

# hglm
fit <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="off")
coef1 <- c(3.04452, -0.45426, -0.29299, 0, 0)
expect_equal(unname(coef(fit)), coef1, tolerance=1e-4)

# predict
model <- hglm(counts~., family=poisson(), data=dobson, constraints=NULL, scaler="off")
gmodel <- glm(counts~., family=poisson(), data=dobson)
expect_equal(predict(model, type="link"), predict(gmodel, type="link"), tolerance=1e-4)
expect_equal(predict(model, newdata=dobson, type="link"), predict(gmodel, newdata=dobson, type="link"), tolerance=1e-4)
expect_equal(predict(model, type="response"), predict(gmodel, type="response"), tolerance=1e-4)
expect_equal(predict(model, newdata=dobson, type="response"), predict(gmodel, newdata=dobson, type="response"), tolerance=1e-4)

# scale predictors
fit_center_standardization <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="center_standardization")
fit_center_minmax <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="center_minmax")
fit_standardization <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="standardization")
fit_minmax <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="minmax")
fit_off <- hglm(counts~., family=poisson(), dobson, constraints=NULL, scaler="off")
coef1 <- c(3.04452, -0.45426, -0.29299, 0, 0)
expect_equal(unname(coef(fit_center_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(fit_center_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(fit_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(fit_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(fit_off)), coef1, tolerance=1e-4)

# scale response + predictors
fit_log_y <- hglm(counts~., family=poisson(link="log"), data=dobson, scale_response=TRUE)
fit_id_y <- hglm(counts~., family=poisson(link="identity"), data=dobson, scale_response=TRUE)
fit_sqrt_y <- hglm(counts~., family=poisson(link="sqrt"), data=dobson, scale_response=TRUE)
expect_equal(unname(coef(fit_log_y)), coef1, tolerance=1e-4)
expect_equal(unname(coef(fit_id_y)), c(21.5307, -7.76271, -5.38842, -0.59051, -0.85045), tolerance=1e-4)
expect_equal(unname(coef(fit_sqrt_y)), c(4.61421, -0.93424, -0.62635, -0.03605, -0.05436), tolerance=1e-4)

## no intercept
# scale predictors
coef_noI <- c(3.04452243772342, 2.59026716544583, 2.75153531304195, 0, 0)
expect_warning(fit_noI_center_standardization <- hglm(counts~.-1, family=poisson(), dobson, constraints=NULL, scaler="center_standardization"), pattern="Intercept .* deactivated")
expect_warning(fit_noI_center_minmax <- hglm(counts~.-1, family=poisson(), dobson, constraints=NULL, scaler="center_minmax"), pattern="Intercept .* deactivated")
fit_noI_standardization <- hglm(counts~.-1, family=poisson(), dobson, constraints=NULL, scaler="standardization")
fit_noI_minmax <- hglm(counts~.-1, family=poisson(), dobson, constraints=NULL, scaler="minmax")
fit_noI_off <- hglm(counts~.-1, family=poisson(), dobson, constraints=NULL, scaler="off")
expect_equal(unname(coef(fit_noI_center_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_center_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_off)), coef_noI, tolerance=1e-4)

## poisson sqrt link
fo <- hp ~ -1 + drat + cyl
coef_noI <- c(0.750701803364248, 1.48388957383876)
expect_warning(fit_noI_center_standardization <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scaler="center_standardization"), pattern="Intercept .* deactivated")
expect_warning(fit_noI_center_minmax <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scaler="center_minmax"), pattern="Intercept .* deactivated")
fit_noI_standardization <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scaler="standardization")
fit_noI_minmax <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scaler="minmax")
fit_noI_off <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scaler="off")
fit_noI_scaleY <- hglm(fo, family=poisson("sqrt"), data=mtcars, constraints=NULL, scale_response=TRUE, scaler="standardization")
expect_equal(unname(coef(fit_noI_center_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_center_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_off)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(fit_noI_scaleY)), coef_noI, tolerance=1e-4)

suppressPackageStartupMessages(library("holiglm"))
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
suppressPackageStartupMessages(library("ROI"))
applicable_solvers <- c("ecos", "gurobi", "cplex")
installed_solvers <- ROI_installed_solvers()
solver <- sample(installed_solvers[names(installed_solvers) %in% applicable_solvers], 1L)
suppressMessages(require(solver, character.only = TRUE))

## Likelihood gaussian

## gaussian identy link function
set.seed(0)
data <- rhglm(100, 1:3, as_list = TRUE)
coef1 <- c(1.06452724839333, 2.01922867737706, 3.08020596165091)
coef_noI <- c(2.07613553650024, 3.11928927313356)

### link function
model <- holiglm:::loglike_gaussian_identity(data$x, data$y)
sol <- ROI_solve(model)
coef0 <- unname(head(solution(sol), NCOL(data$x)))
expect_equal(coef0, coef1, tolerance=1e-4)

### hglm_model
model <- hglm_model(data$x, data$y, family=gaussian())
op <- as.OP(model)
sol <- ROI_solve(op)
coef0 <- unname(head(solution(sol), NCOL(data$x)))
expect_equal(coef0, coef1, tolerance=1e-4)

### hglm
df <- as.data.frame(data$x[,-1])
df$y <- data$y
model <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="off")
coef0 <- unname(coef(model))
expect_equal(coef0, coef1, tolerance=1e-4)

# scaling predictors
model_center_standardization <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="center_standardization")
model_center_minmax <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="center_minmax")
model_standardization <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="standardization")
model_minmax <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="minmax")
model_off <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="off")
expect_equal(unname(coef(model_center_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_center_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_off)), coef1, tolerance=1e-4)

# scaling response + predictors
model_y_center_standardization <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="center_standardization", scale_response=TRUE)
model_y_center_minmax <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="center_minmax", scale_response=TRUE)
model_y_standardization <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="standardization", scale_response=TRUE)
model_y_minmax <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="minmax", scale_response=TRUE)
model_y_off <- hglm(y~., family=gaussian(), df, constraints=NULL, scaler="off", scale_response=TRUE)
expect_equal(unname(coef(model_y_center_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_y_center_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_y_standardization)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_y_minmax)), coef1, tolerance=1e-4)
expect_equal(unname(coef(model_y_off)), coef1, tolerance=1e-4)

## no intercept
# scaling predictors
expect_warning(model_noI_center_standardization <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="center_standardization"), pattern="Intercept .* deactivated")
expect_warning(model_noI_center_minmax <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="center_minmax"), pattern="Intercept .* deactivated")
model_noI_standardization <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="standardization")
model_noI_minmax <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="minmax")
model_noI_off <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="off")
expect_equal(unname(coef(model_noI_center_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_noI_center_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_noI_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_noI_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_noI_off)), coef_noI, tolerance=1e-4)

# scaling response + predictors
expect_warning(model_y_noI_center_standardization <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="center_standardization", scale_response=TRUE), pattern="Intercept .* deactivated")
expect_warning(model_y_noI_center_minmax <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="center_minmax", scale_response=TRUE), pattern="Intercept .* deactivated")
model_y_noI_standardization <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="standardization", scale_response=TRUE)
model_y_noI_minmax <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="minmax", scale_response=TRUE)
model_y_noI_off <- hglm(y~.-1, family=gaussian(), df, constraints=NULL, scaler="off", scale_response=TRUE)
expect_equal(unname(coef(model_y_noI_center_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_y_noI_center_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_y_noI_standardization)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_y_noI_minmax)), coef_noI, tolerance=1e-4)
expect_equal(unname(coef(model_y_noI_off)), coef_noI, tolerance=1e-4)

## gaussian log link function
# set.seed(0)
# true_beta <-c(0.8, 0.3)
# data <- rhglm(1000, true_beta, as_list = TRUE, family = gaussian("log"))
# model <- holiglm:::loglike_gaussian_log(data$x, data$y)
# sol <- ROI_solve(model)
# coef0 <- NULL
# coef1 <- unname(head(solution(sol, force = TRUE), NCOL(data$x)))
# coef2 <- unname(coef(glm.fit(data$x, data$y, family = gaussian("log"), start = true_beta)))
# rbind(coef1, coef2)
# cat(deparse(coef2))
# expect_equal(coef0, coef1, tolerance=1e-4)

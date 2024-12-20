suppressPackageStartupMessages(library("holiglm"))
## Check auxiliary functions

## Check if link is power function
expect_true( holiglm:::is_power(binomial(link=power(2))))
expect_true( holiglm:::is_power(poisson(link=power(1/2))))
expect_false(holiglm:::is_power(binomial()))
expect_false(holiglm:::is_power(gaussian()))
expect_false(holiglm:::is_power(poisson()))
# link=power(0) changes link to link="log"
expect_false(holiglm:::is_power(gaussian(link=power(0))))
# link=power(1) changes link to link="identity"
expect_false(holiglm:::is_power(gaussian(link=power(1))))


## Extract power from power link
expect_equal(holiglm:::extract_power(poisson(link=power(2))), 2)
expect_equal(holiglm:::extract_power(poisson(link=power(1/2))), round(1/2, 3))
expect_equal(holiglm:::extract_power(binomial(link=power(1/3))), round(1/3, 3))
expect_error(holiglm:::extract_power(gaussian()), "Assertion on.*family.*")
expect_error(holiglm:::extract_power(binomial()), "Assertion on.*family.*")
expect_error(holiglm:::extract_power(poisson()), "Assertion on.*family.*")


## Convert two level factor to binary variable
f1 <- iris$Species[iris$Species != "setosa"]
expect_equal(holiglm:::factor_binary(f1), as.integer(f1)-2L)

# factor contains more than two levels
expect_error(holiglm:::factor_binary(iris$Species))
# not a factor
expect_error(holiglm:::factor_binary(1:10))

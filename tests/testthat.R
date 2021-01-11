library(testthat)
library(MiRKATMC)

test_check("MiRKATMC")

test.data <- data.frame(y = rep(c(1:4), 20), x1 = rnorm(80, 0, 1), x2 = rbinom(80, 1, 0.5))

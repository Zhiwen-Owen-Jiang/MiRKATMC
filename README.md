# MiRKATMC
The Microbiome Regression-based Kernel Association Test for multi-categorical (nominal or ordinal) data

## Installation
```{r}
# install.packages("devtools") if you have not installed devtools
devtools::install_github("Zhiwen-Owen-Jiang/MiRKATMC")
```
## Usage
```{r}
library(MiRKATMC)

# let's first generate some data
set.seed(123)
test.data <- data.frame(outcome = as.factor(sample(4, 100, replace = TRUE)),
                        ID = gl(20, 5), time = rep(1:5, 20), age = rnorm(n = 100, mean = 30, sd = 5),
                        sex = rbinom(100, 1, 1/2))
D1 <- matrix(rbinom(10000, 2, 0.05), 100, 100)
K1 <- crossprod(D1) # kernel matrix
D2 <- matrix(rbinom(10000, 2, 0.1), 100, 100)
K2 <- crossprod(D2) # kernel matrix
K <- list(kernel1 = K1, kernel2 = K2)
K_no_name <- list(K1, K2)

# Then do the analysis
MiRKATMC(formula = outcome ~ age, random = NULL, data.type = 'nominal', Ks = K1, data = test.data)
MiRKATMC(formula = outcome ~ age, random = NULL, data.type = 'ordinal', Ks = K_no_name, data = test.data)
MiRKATMC(formula = outcome ~ age, random = ~ 1 | ID, data.type = 'nominal', Ks = K, data = test.data)
MiRKATMC(formula = outcome ~ age, random = ~ 1 + time | ID, data.type = 'ordinal', Ks = K_no_name, data = test.data)
```
## Getting help
Please email Owen Jiang <owenjf@live.unc.edu>.

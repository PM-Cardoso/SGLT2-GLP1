####################
## Description:
##  - In this file we make a simulation study.
##    Goal: Investigate whether grf and bcf can calculate treatment effect well
##          as well as investigate if bartMachine can do it too without being causal.
##
##  - How:
##        - Model1: - X1:3 from a multivariate Gaussian
##                  - Z: binary treatment effect
##                  - y: Z + Z*X + error
##        - Model2: - X1:6 from a multivariate Gaussian
##                  - Z: binary treatment effect from X1:3 (confounding)
##                  - y: Z + Z*X + error
####################

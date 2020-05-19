rm(list=ls(all=TRUE))

mu_true = 10
n = 100
a0 = 0.01
g0 = 0.01
a = a0 + n
g = g0 + sum(y-mu_true)^2
h1 = rgamma(1,a/2,g/2)

library(glmnet)
library(ISLR)

data(Hitters, package = "ISLR")
Hitters = na.omit(Hitters)


a = c(1,2,3)
b = c(2,4,1)

library(Matrix)
library(geometry)
library(dplyr)


m1 <- Matrix(nrow = 1, ncol = 1000, data = 0, sparse = TRUE)
m1 <- as(m, "dgTMatrix") # by default, Matrix() returns dgCMatrix
m1


m2 <- Matrix(nrow = 1, ncol = 1000, data = 0, sparse = TRUE)
m2 <- as(m, "dgTMatrix") # by default, Matrix() returns dgCMatrix
m2


matSparse <- as(matrix(data = c(1,0,0,0,3,0), nrow = 1), "dgCMatrix")
matSparse

scalar = as(matrix(data = c(1,1,1,1,1,0), nrow = 1), "dgCMatrix")
scalar

matSparse * scalar %>% as_tibble()

matSparse %*% t(scalar)

geometry::dot(matSparse , scalar)



t1<-Sys.time()

sv1 <- sparseVector(x = c(1,2,3), i = c(1,3,5), length = 1000)
sv1

sv2 <- sparseVector(x = c(1,1,1), i = c(1,3,5), length = 1000)
sv2

sv3 = sv2*t(sv1) 
#sv3@x %>% sum
sv3@p[length(sv3@p)]

t2<-Sys.time()
print(t2 - t1)
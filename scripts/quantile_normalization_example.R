# https://www.bioconductor.org/help/course-materials/2019/NURDS/NURDS.html
x <- runif(5)
x

xo <- x[order(x)]
xo

xo[order(order(x))]

diff(c(0, xo))[order(order(x))]

LETTERS[seq_along(xo)][order(order(x))]

set.seed(123)
n_genes <- 1000
base <- 2 + rexp(n_genes)
base

expression <- data.frame(
  x = base + rnorm(n_genes, 0.3) / base,
  y = 1.1 * base + rnorm(n_genes, 0.4) / base
)

library(ggplot2)
ggplot(expression, aes(x, y)) + geom_point()

df <- reshape2::melt(expression)

ggplot(df, aes(x = value, colour = factor(variable))) +
  geom_density()

x <- expression$x
xo <- x[order(x)]
y <- expression$y
yo <- y[order(y)]

mo <- cbind(xo, yo)
row_mean <- apply(mo, 1, mean)

normalized <- data.frame(
  x = row_mean[order(order(x))],
  y = row_mean[order(order(y))]
)

normalized

ggplot(normalized, aes(x, y)) + geom_point()

df <- reshape2::melt(normalized)

ggplot(df, aes(x = value, colour = factor(variable))) +
  geom_density()

quantile_normalize <-
  function(expression) {
    x <- expression$x
    xo <- x[order(x)]
    y <- expression$y
    yo <- y[order(y)]
    
    mo <- cbind(xo, yo)
    row_mean <- apply(mo, 1, mean)
    
    data.frame(
      x = row_mean[order(order(x))],
      y = row_mean[order(order(y))]
    )
  }

normalized <- quantile_normalize(expression)

expression <- data.frame(
  x = c(1, 2, 4, 7),
  y = c(3, 2, 5, 8)
)
quantile_normalize(expression)

library(testthat)
test_that("quantile_normalize() works", {
  ## calculate simple example 'by hand'...
  x <- c(1, 2, 4, 7)
  y <- c(3, 2, 5, 8)
  xo <- x[order(x)]
  yo <- y[order(y)]
  mo <- apply(cbind(xo, yo), 1, mean)
  expected <- data.frame(x = mo[order(order(x))], y = mo[order(order(y))])
  
  ## Compare to outcome of our function
  expression <- data.frame(x = x, y = y)
  observed <- quantile_normalize(expression)
  expect_equal(observed, expected)      # other expectations possible
})

quantile_normalize <-
  function(x) {
    m <- apply(x, 2, function(v) v[order(v)])
    row_mean <- apply(m, 1, mean)
    apply(x, 2, function(v, row_mean) row_mean[order(order(v))], row_mean)
  }

quantile_normalize <-
  function(x){
    m <- apply(x, 2, function(v) v[order(v)])
    row_mean <- rowMeans(m)
    apply(x, 2, function(v, row_mean) row_mean[order(order(v))], row_mean)
  }

test_that("quantile_normalize() works", {
  x <- c(1, 2, 4, 7)
  y <- c(3, 2, 5, 8)
  xo <- x[order(x)]
  yo <- y[order(y)]
  mo <- apply(cbind(xo, yo), 1, mean)
  expected <- cbind(x = mo[order(order(x))], y = mo[order(order(y))])
  
  expression <- data.frame(x = x, y = y)
  observed <- quantile_normalize(expression)
  expect_equal(observed, expected)
})

quantile_normalize <-
  function(x) {
    ## so long as the input can be coerced to a matrix...
    x <- as.matrix(x)
    
    ## and can be validated to conform to our contract...
    stopifnot(
      is.numeric(x),
      !anyNA(x)
    )
    
    ## ...we have confidence that we will satisfy the return value
    m <- apply(x, 2, function(v) v[order(v)])
    row_mean <- rowMeans(m)
    apply(x, 2, function(v, mo) mo[order(order(v))], mo)
  }

test_that("'quantile_normalize()' validates inputs", {
  m <- cbind(letters, LETTERS)
  expect_error(quantile_normalize(m))
  
  df <- data.frame(x=rnorm(26), y = letters)
  expect_error(quantile_normalize(df))
  
  m <- matrix(rnorm(10), nrow = 5)
  m[1,1] <- NA
  expect_error(quantile_normalize(m))
})

quantile_normalize <-
  function(x) {
    ## validate inputs
    x <- as.matrix(x)
    stopifnot(
      is.numeric(x),
      !anyNA(x)
    )
    
    ## quantile normalization
    m <- apply(x, 2, function(v) v[order(v)])
    row_mean <- rowMeans(m)
    result <-
      apply(x, 2, function(v, row_mean) row_mean[order(order(v))], row_mean)
    
    ## propagate dimnames
    dimnames(result) <- dimnames(x)
    result
  }

test_that("'quantile_normalize()' propagates dimnames", {
  m <- matrix(rnorm(10), 5, dimnames=list(LETTERS[1:5], letters[1:2]))
  observed <- quantile_normalize(m)
  expect_identical(dimnames(observed), dimnames(m))
})

quantile_normalize <-
  function(x) {
    ## validate inputs
    x <- as.matrix(x)
    stopifnot(
      is.numeric(x),
      !anyNA(x)
    )
    
    ## quantile normalize
    m <- apply(x, 2, function(v) v[order(v)])
    dim(m) <- dim(x)    # apply() doesn't always return a matrix!
    
    row_mean <- rowMeans(m)
    
    result <- apply(
      x, 2, function(v, row_mean) row_mean[order(order(v))], row_mean
    )
    dim(result) <- dim(x)
    
    ## propagate dimnames
    dimnames(result) <- dimnames(x)
    result
  }

test_that("'quantile_normalize()' works with edge cases", {
  m <- matrix(rnorm(5), nrow = 5)
  expect_identical(m, quantile_normalize(m))
  
  m <- matrix(rnorm(5), ncol = 5)
  expect_identical(matrix(mean(m), ncol = 5), quantile_normalize(m))
  
  m <- matrix(0, 0, 0)
  expect_identical(m, quantile_normalize(m))
})

library(microbenchmark)

n_genes <- 10000
g10 <- matrix(rnorm(n_genes * 10), ncol = 10)
g100 <- matrix(rnorm(n_genes * 100), ncol = 100)
g1000 <- matrix(rnorm(n_genes * 1000), ncol = 1000)
g10000 <- matrix(rnorm(n_genes * 10000), ncol = 10000)

times <- microbenchmark(
  quantile_normalize(g10),
  quantile_normalize(g100),
  quantile_normalize(g1000),
  quantile_normalize(g10000),
  times = 10
)

times

Rprof()        # start profiling
result <- quantile_normalize(g10000)
Rprof(NULL)    # stop profiling
summaryRprof()$by.total

# some issues in multi-threading, reinstall the package with disabling it
# https://stackoverflow.com/questions/61629861/error-return-code-from-pthread-create-is-22
# BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
library(preprocessCore)
normalize.quantiles(g10000)

microbenchmark(quantile_normalize(g10), normalize.quantiles(g10), times = 10)



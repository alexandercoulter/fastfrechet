test_that("generate_zinbinom_qf output meets expected structure and properties", {
  # Parameters
  n <- 100  # Number of samples (rows of X and Y)
  p <- 10   # Number of covariates (columns of X)
  m <- 100  # grid density (columns of Y)
  
  # Generate data
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X  # (n x p) matrix of covariates
  Y <- mydata$Y  # (n x m) matrix of empirical quantile functions, stored row-wise
  
  # Explicitly check rows and columns of X
  expect_equal(nrow(X), n, info = "Number of rows in X is incorrect.")
  expect_equal(ncol(X), p, info = "Number of columns in X is incorrect.")
  
  # Explicitly check rows and columns of Y
  expect_equal(nrow(Y), n, info = "Number of rows in Y is incorrect.")
  expect_equal(ncol(Y), m, info = "Number of columns in Y is incorrect.")
  
  # Check that each row of Y is non-decreasing
  non_decreasing_check <- apply(Y, 1, function(row) all(diff(row) >= 0))
  expect_true(all(non_decreasing_check), 
              info = "Not all rows of Y are non-decreasing, indicating invalid quantile functions.")
})
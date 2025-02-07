test_that("predicted QFs are monotonic and within [0, 40]", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = NULL,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with large p", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = X,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with lambda", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Generate random lambda vector with \tau sum
  tau <- 10
  g <- rnorm(p)
  g <- g / sqrt(sum(g * g))
  lambda <- g * g * tau
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = X,
                                   C_init = NULL,
                                   lambda = lambda,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with lambda and large p", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Generate random lambda vector with \tau sum
  tau <- 10
  g <- rnorm(p)
  g <- g / sqrt(sum(g * g))
  lambda <- g * g * tau
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = X,
                                   C_init = NULL,
                                   lambda = lambda,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with nonsense C_init", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Estimate conditional QFs with random initial active constraint matrix:
  output1 <- frechetreg_univar2wass(X = X,
                                    Y = Y,
                                    Z = X,
                                    C_init = matrix(rbinom(n * (m + 1), 1, 0.5), n, m + 1),
                                    lambda = NULL,
                                    lower = lower,
                                    upper = upper)
  # Estimate conditional QFs with default active constraint matrix (all zeros):
  output2 <- frechetreg_univar2wass(X = X,
                                    Y = Y,
                                    Z = X,
                                    C_init = NULL,
                                    lambda = NULL,
                                    lower = lower,
                                    upper = upper)
  
  # Extract predicted QFs and take difference:
  predicted_QF_diff <- output1$Qhat - output2$Qhat # Assuming these contain (n x m) matrices of QFs
  
  # Check L-infinity norm is bounded to high tolerance:
  expect_true(max(abs(predicted_QF_diff)) <= 1e-10, 
              info = "QF values violate the lower bound")
})


test_that("works with different output matrix Z", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Generate Z matrix
  nz = 10
  Z <- matrix(rnorm(nz * p), nz, p)
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = Z,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with different output matrix Z, and large p", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Generate Z matrix
  nz = 10
  Z <- matrix(rnorm(nz * p), nz, p)
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = Z,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with different output matrix Z, large p, and lambda", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] <- upper
  
  # Generate Z matrix
  nz = 10
  Z <- matrix(rnorm(nz * p), nz, p)
  
  # Generate random lambda vector with \tau sum
  tau <- 10
  g <- rnorm(p)
  g <- g / sqrt(sum(g * g))
  lambda <- g * g * tau
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = Z,
                                   C_init = NULL,
                                   lambda = lambda,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})

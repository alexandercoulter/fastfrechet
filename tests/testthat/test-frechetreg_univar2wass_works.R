test_that("predicted QFs are monotonic and within [0, 40]", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with large p", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with large n but non full column rank X", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- cbind(mydata$X, mydata$X)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with Z", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  nz <- 10
  Z <- matrix(rnorm(nz * p), nz, p)

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("fails with large p and improper Z", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  nz <- 10
  Z <- matrix(rnorm(nz * p), nz, p)

  expect_error(frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  ), "The rows of matrix (1 Z) are not all in the row-space of (1 X).", fixed = TRUE)
})


test_that("works with large p and Z passing row-space condition", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  nz <- 10
  A <- matrix(rnorm(nz * n), nz, n)
  A <- A / rowSums(A)
  Z <- A %*% X

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("fails with arbitrary Z, and large n but non full column rank X", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- cbind(mydata$X, mydata$X)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Generate Z matrix
  nz <- 10
  Z <- matrix(rnorm(nz * p * 2), nz, p * 2)

  expect_error(frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  ), "The rows of matrix (1 Z) are not all in the row-space of (1 X).", fixed = TRUE)
})


test_that("works with Z satisfying row-space condition, and large n but non full column rank X", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- cbind(mydata$X, mydata$X)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Generate Z matrix
  nz <- 10
  A <- matrix(rnorm(nz * n), nz, n)
  A <- A / rowSums(A)
  Z <- A %*% X

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with lambda", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with lambda and large p", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with lambda and Z", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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

  # Generate Z matrix
  nz <- 10
  Z <- matrix(rnorm(nz * p), nz, p)

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with lambda, Z, and large p", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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

  # Generate Z matrix
  nz <- 10
  Z <- matrix(rnorm(nz * p), nz, p)

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = Z,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with nonsense C_init", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

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
  output1 <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = X,
    C_init = matrix(rbinom(n * (m + 1), 1, 0.5), n, m + 1),
    lambda = NULL,
    lower = lower,
    upper = upper
  )
  # Estimate conditional QFs with default active constraint matrix (all zeros):
  output2 <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = X,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs and take difference:
  predicted_QF_diff <- output1$Qhat - output2$Qhat # Assuming these contain (n x m) matrices of QFs

  # Check L-infinity norm is bounded to high tolerance:
  expect_true(max(abs(predicted_QF_diff)) <= 1e-10,
    info = "QF values differ by C_init"
  )
})


test_that("works with X = 0", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- matrix(0, n, p)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with X = 0 and large p", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- matrix(0, n, p)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = NULL,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with X = 0 and lambda", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- matrix(0, n, p)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Generate random lambda vector with \tau sum
  tau <- 10
  g <- rnorm(p)
  g <- g / sqrt(sum(g * g))
  lambda <- g * g * tau

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})


test_that("works with X = 0, lambda, and large p", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- matrix(0, n, p)
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Generate random lambda vector with \tau sum
  tau <- 10
  g <- rnorm(p)
  g <- g / sqrt(sum(g * g))
  lambda <- g * g * tau

  # Estimate conditional QFs
  output <- frechetreg_univar2wass(
    X = X,
    Y = Y,
    Z = NULL,
    C_init = NULL,
    lambda = lambda,
    lower = lower,
    upper = upper
  )

  # Extract predicted QFs
  predicted_QFs <- output$Qhat # Assuming this contains (n x m) matrix of QFs

  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })

  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5),
    info = "QF values violate the lower bound"
  )

  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5),
    info = "QF values violate the upper bound of 40"
  )
})

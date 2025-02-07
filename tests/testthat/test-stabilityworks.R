# Parameters for the example
n <- 100 # number of samples
p <- 10 # number of covariates
m <- 50 # EQF grid density
lower <- 0
upper <- Inf
B <- 20
thresh <- 0.0001
tauseq <- 1:10
eps <- 0.001

# Generate data
set.seed(31)
mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
X <- mydata$X
Y <- mydata$Y

# Run complementary pairs stability selection
cpss <- FRiSO_CPSS_univar2wass(
  X = X,
  Y = Y,
  B = B,
  thresh = thresh,
  lower = lower,
  upper = upper,
  tauseq = tauseq,
  eps = eps
)



test_that("stability_paths rows are non-increasing", {
  # Extract stability paths
  stability_paths <- cpss$stability_paths

  apply(stability_paths, 2, function(column) {
    # Check that last row value > first row value
    expect_true(column[length(column)] >= column[1],
      info = "Largest tau gives smaller stability probability than smallest tau"
    )
    # Check that all values are strictly between 0 and 1
    expect_true(all(column >= 0 & column <= 1),
      info = "Stability probabilities outside of [0,1]"
    )
  })
})

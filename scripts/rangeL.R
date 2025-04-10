# Authenticate and install package if not installed
if (!requireNamespace("fastfrechet", quietly = TRUE)) {
  remotes::install_github("alexandercoulter/fastfrechet", auth_token = Sys.getenv("GITHUB_PAT"))
}
library(fastfrechet)

n <- 200 # number of samples - nrow(X) and nrow(Y).
p <- 20 # number of covariates - ncol(X).
m <- 200 # quantile function grid density - ncol(Y).
mseq <- seq(1 / (2 * m), 1 - 1 / (2 * m), length.out = m)

set.seed(31)
mydata <- generate_zinbinom_qf(
  n = n,
  p = p,
  m = m
)

X <- mydata$X # (n x p) matrix of covariates
Y <- mydata$Y # (n x m) matrix of quantile functions, stored row-wise

# Dense grid of "allowance" totals:
tauseq <- seq(0.2, 20, 0.2)

# Generate estimated "allowance vector"s \lambda for each \tau, stored
# column-wise in matrix `L`:
L <- FRiSO_univar2wass(
  X = X,
  Y = Y,
  lower = 0,
  upper = Inf,
  tauseq = tauseq,
  eps = 0.001,
  nudge = 0.01
)

# Calculate range of L:
rangeL = range(L)

# Create output directory:
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)

write.csv(L, file = file.path(output_dir, "L.csv"))
write.csv(rangeL, file = file.path(output_dir, "rangeL.csv"))

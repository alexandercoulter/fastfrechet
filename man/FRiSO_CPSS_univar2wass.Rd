% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FRiSO_CPSS_univar2wass.R
\name{FRiSO_CPSS_univar2wass}
\alias{FRiSO_CPSS_univar2wass}
\title{Complementary pairs stability selection (CPSS) for FRiSO, for univariate distributions with 2-Wasserstein metric}
\usage{
FRiSO_CPSS_univar2wass(X, Y, B = 50, thresh = 1e-04, ...)
}
\arguments{
\item{X}{A (\code{n} \eqn{\times} \code{p}) "input" covariate matrix with no missing, all finite entries. We recommend you \strong{do not} include an intercept column (i.e. all \code{1}'s vector).}

\item{Y}{A (\code{n} \eqn{\times} \code{m}) matrix of observed quantile functions evaluated on a shared uniform \code{m}-grid in \eqn{(0, 1)} - e.g. \code{seq(0.5/m, 1 - 0.5/m, len = m)} - row-wise monotone non-decreasing. Entries must obey user-specified box constraints given by \code{lower} and \code{upper} parameters.}

\item{B}{A positive integer number of complementary pair sub-samples to generate from the original data set (default \code{50}).}

\item{thresh}{A finite positive scalar selection threshold, where \eqn{\widehat{\lambda}_k(\tau) > } \code{thresh} means the \eqn{k^{\mathrm{th}}} variable is selected; otherwise it is not selected.}

\item{...}{other inputs to \code{\link[=FRiSO_univar2wass]{FRiSO_univar2wass()}}.}
}
\value{
\tabular{ll}{
\code{tau} \tab returns a vector containing the values of \eqn{\tau} which were fed into the function call. \cr
\code{selected_variables} \tab returns a (\code{B} \eqn{\times} \code{2} \eqn{\times} \code{length(tauseq)} \eqn{\times} \code{n}) array containing 0's and 1's identifying which samples were selected in course of CPSS.\cr
\code{selected_samples} \tab returns a (\code{B} \eqn{\times} \code{2} \eqn{\times} \code{length(tauseq)} \eqn{\times} \code{n}) array containing 0's and 1's identifying which samples were selected in course of CPSS. \cr
\code{stability_paths} \tab returns a (\code{length(tauseq)} \eqn{\times} \code{p}) matrix containing stability measures for each variable (column-wise) against given tauseq.
}
}
\description{
This function performs complementary pairs (\insertCite{shah_variable_2013}{fastfrechet})
stability selection (\insertCite{meinshausen_stability_2010}{fastfrechet})
with the Fréchet Ridge Selection Operator (FRiSO; \insertCite{tucker_variable_2023}{fastfrechet}),
for the space of univariate distribution responses equipped with the 2-Wasserstein
metric. Heuristically, the optimal model is identified by selecting covariates
whose "stability paths" over \eqn{\tau} exceed a threshold, which may itself
be a function of \eqn{\tau}.
}
\examples{
# Generate data for X and Y inputs by using the output of `generate_zinbinom_qf`
n <- 100 # number of samples - nrow(X) and nrow(Y).
p <- 10 # number of covariates - ncol(X).
m <- 50 # EQF grid density - ncol(Y).
lower <- 0
upper <- Inf

set.seed(31)
mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
X <- mydata$X # (n x p) matrix of covariates
Y <- mydata$Y # (n x m) matrix of EQFs, stored row-wise

# Set complementary pairs stability selection parameters
B <- 50
thresh <- 0.0001
tauseq <- 1:10
eps <- 0.001

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

# Plot stability paths
matplot(cpss$tau, cpss$stability_paths,
  type = "l", lty = 1, lwd = 2,
  col = c(rep("red", 4), rep("black", p - 4))
)

# Calculate thresholds using Shah and Samworth (2013) method:
shahsam <- Shah_Samworth_thresholds(
  p = p,
  q = cpss$model_size_est,
  B = B,
  E_thr = 1
)

# Add lines to plot:
lines(cpss$tau, shahsam$pointwise_thresholds, type = "o", pch = 16)
}
\references{
\insertRef{shah_variable_2013}{fastfrechet}

\insertRef{meinshausen_stability_2010}{fastfrechet}

\insertRef{tucker_variable_2023}{fastfrechet}
}

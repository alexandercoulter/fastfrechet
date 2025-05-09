% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monotoneQP.R
\name{monotoneQP}
\alias{monotoneQP}
\title{Wrapper for QP solver for nearest empirical quantile function in Frobenius norm}
\usage{
monotoneQP(Y, C_init = NULL, lower = -Inf, upper = Inf, eps = 1e-10)
}
\arguments{
\item{Y}{An (\code{m})-long numeric vector (or (\code{n} \eqn{\times} \code{m}) matrix) with no missing, all finite entries.}

\item{C_init}{An optional (\code{m + 1})-long numeric vector (or (\code{n} \eqn{\times} \code{m + 1}) matrix) of non-negative entries, specifying initial active set(s) for optimization. Active sets are identified by positive entries, row-wise if a matrix.}

\item{lower}{An optional numeric scalar (default \code{-Inf}) lower box constraint; must be strictly less than \code{upper}.}

\item{upper}{An optional numeric scalar (default \code{Inf}) upper box constraint; must be strictly greater than \code{lower}.}

\item{eps}{A positive numeric scalar error tolerance (default \code{1e-10}).}
}
\value{
A list object with components:
\tabular{ll}{
\code{Lagrange_Multiplier} \tab returns a (\code{n} \eqn{\times} (\code{m + 1})) matrix, row-wise containing the Lagrange multipliers associated with the QP problems. \cr
\code{Solution} \tab returns a (\code{n} \eqn{\times} \code{m}) matrix, row-wise containing the solutions to the QP problems. \cr
}
}
\description{
\code{Rcpp} implementation of a custom dual active-set algorithm for finding the nearest
discretized quantile function in Frobenius norm to an input vector. The function
permits a matrix input, where the optimization problem is solved across the
rows of the matrix. is solved row-wise, and also permits warm starts, in
the form of specifying an initial active-constraint matrix. Implementation
is based on the dual-active set method of
\insertCite{arnstrom_dual_2022}{fastfrechet}, taking advantage of
simplifications in this setting which avoid \eqn{LDL^{\top}} decompositions
and other costly matrix product operations.
}
\details{
This function solves \deqn{\widehat{\mathbf{q}} := \underset{\mathbf{q} \in \mathbb{R}^m}{\mathrm{argmin}}\lVert \mathbf{q} - \mathbf{y} \rVert_2^2}
for input vector \eqn{\mathbf{y}}, with the constraints \eqn{b_L \leq q_1 \leq \dots \leq q_m \leq b_U} for \eqn{b_L < b_U}.
The fitted vector \eqn{\widehat{\mathbf{q}}} is thus a discretized quantile function
with box constraints, which can be found through quadratic programming. This
function implements a customized dual active-set method inspired by
\insertCite{arnstrom_dual_2022}{fastfrechet}, which takes advantage of
sparsity in the constraints to avoid \eqn{LDL^{\top}} decompositions and other
costly matrix operations.

The user can specify a warm start for the algorithm in the form of an \code{m + 1}-long
vector (or an \code{n} \eqn{\times} \code{m + 1} matrix, if \code{Y} is an \code{n} \eqn{\times} \code{m}
matrix) which gives an estimate of the solution's active set(s), i.e. which
of the \code{m + 1} constraints per row of \code{Y} are exact equalities. Positive
entries in the optional input \code{C_init} indicate active constraints, and zero
entries indicate inactive constraints. (If \code{C_init} contains any negative
entries, the function will stop with an error message.) The default is to
assume no constraints are active.

As this function implements an active-set method, the solution should be
exact to numerical precision. For stability, however, it is advised to keep
the error tolerance parameter \code{eps} at a very small, positive number, such
as the default \code{1e-10}.
}
\examples{
# Set box constraints:
lower <- 0.5
upper <- 1.5

# Generate example vector:
m <- 100
set.seed(31)
y <- rnorm(m, 2 * seq(0, 1, len = m), 0.1)

# Calculate monotone, box-constrained projection:
output <- monotoneQP(y, lower = lower, upper = upper)

# Note: to numerical precision, these QFs are non-decreasing...
min(diff(output$Solution[1, ]))

# ...bounded from below by the lower bound...
min(output$Solution) - lower

# ...and bounded from above by the upper bound...
max(output$Solution) - upper

# Plot values of the generated vector:
plot(y,
  main = "Monotone and Box Constrained", las = 1,
  xlab = "Vector Entry", ylab = expression("y, " * hat("q"))
)
abline(h = c(lower, upper), lty = 2, col = "gray80")

# Add values of monotone projection:
points(output$Solution[1, ], pch = 20)

legend("topleft",
  pch = c(1, 20), bty = "n",
  legend = c("Unconstrained", "Monotone Constrained")
)
}
\references{
\insertRef{arnstrom_dual_2022}{fastfrechet}
}

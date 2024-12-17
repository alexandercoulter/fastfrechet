---
title: 'fastfrechet: An R package for fast implementation of Fréchet regression functions
  for distribution responses'
tags:
- Fréchet regression
- Variable selection
- Distributions
- Wasserstein
- Metric space
- R
- C++
date: "27 November 2024"
affiliations:
- name: Department of Statistics, Texas A&M University, United States
  index: 1
- name: Department of Biostatistics, University of Michigan, United States
  index: 2
authors:
- name: Alexander Coulter
  orcid: "0009-0001-9116-6588"
  equal-contrib: yes
  affiliation: 1
- name: Rebecca Lee
  orcid: "0000-0001-8162-6555"
  equal-contrib: yes
  affiliation: 1
- name: Dr. Irina Gaynanova
  orcid: "0000-0002-4116-0268"
  equal-contrib: yes
  affiliation: 2
  corresponding: yes
bibliography: inst/REFERENCES.bib
---

# Summary

Distribution-as-response regression problems are gaining wider attention,
especially within biomedical settings where observation-rich patient specific
data sets are available, such as continuous glucose monitoring
[@matabuena_glucodensities_2021], actigraphy [@ghosal_distributional_2023], and
feature densities in CT scans [@petersen_wasserstein_2021].
To accommodate the complex structure of such problems, @petersen_frechet_2019
proposed an extension of regression with Euclidean covariates, for responses
within a general metric space, so-called *Fréchet regression*. This regression
framework was further extended for variable selection by @tucker_variable_2023,
and @coulter_fast_2024 developed a fast variable selection algorithm for the
specific setting of univariate distribution responses equipped with the
2-Wasserstein metric (*2-Wasserstein space*). We present `fastfrechet`, an R
package providing fast implementation of these Fréchet regression and variable
selection methods in 2-Wasserstein space, with resampling tools for automatic
variable selection. `fastfrechet` makes 2-Wasserstein Fréchet regression with
resampling-supplemented variable selection readily available and highly scalable
to large data sets, such as the UK Biobank [@doherty_large_2017].

# Statement of Need

No software package currently supports variable selection for 2-Wasserstein
Fréchet regression, and until @coulter_fast_2024, existing algorithms are too
slow for practical use. To our knowledge, two packages implement Fréchet
regression in 2-Wasserstein space: the R package `WRI`, and the R package
`frechet`. However, these packages are computationally inefficient and limited
in scope. For example, `WRI` does not permit user-specified box constraints
for distribution support, and requires strictly increasing quantile function
inputs. `frechet` allows more flexible user specifications, but its Fréchet
regression solver is slow and does not satisfy constraints to satisfactory
accuracy.

`fastfrechet` bridges these gaps by providing a fast and scalable package-ready
implementation of variable selection for 2-Wasserstein Fréchet regression. It
includes resampling tools—cross-validation as discussed in
@tucker_variable_2023, and stability selection as discussed in
@coulter_fast_2024—to supplement automatic variable selection. `fastfrechet`
also implements a customized dual active-set solver based on @arnstrom_dual_2022
to solve the Fréchet regression problem, which overcomes computational and user
specification limitations in existing packages.

# Performance Comparisons to Existing Implementations

We illustrate the performance of R package `fastfrechet` against existing
implementations with simulated covariate-dependent, zero-inflated negative
binomial distributions. The included function `generate_zinbinom_qf` allows us
to simulate $n$ such distributions, represented as quantile functions evaluated
on a shared $m$-grid in $(0, 1)$, dependent on the first 4 of $p \geq 4$
covariates. For a general illustration of what Fréchet regression and variable
selection look like in 2-Wasserstein space, see the accompanying
`intro-fastfrechet` vignette.
```
set.seed(31)
n = 100          # number of QFs
p = 10           # number of covariates
m = 100          # (0, 1) QF grid density 

gendata = fastfrechet::generate_zinbinom_qf(n, p, m)
X = gendata$X    # (n x p) covariate matrix
Y = gendata$Y    # (n x m) QF response matrix, stored row-wise
```

### The Fréchet Regression Problem

The R package `fastfrechet` provides a solver for the Fréchet regression problem
for 2-Wasserstein space, with optional user-specified box constraints to enforce
possibly finite support on the underlying distributions. The implementation is
a customization of the dual active-set method of @arnstrom_dual_2022 (see the
accompanying `monotoneQP-fastfrechet` vignette). The resulting algorithm
provides a fast and numerically exact solution to the Fréchet regression
problem, illustrated in Figure \autoref{fig:frechetreg_comparison}. Since the
`WRI` package requires strictly monotone inputs, we add a small adjustment to
the response matrix when using `WRI`'s solver.
```
library(microbenchmark)

X_df = as.data.frame(X)
Y_adj = Y + matrix(seq(0, 0.01, len = m), n, m, byrow = TRUE)

microbenchmark(WRI::wass_regress( ~ ., X_df, "quantile", Y_adj),
               frechet::GloDenReg(X, qin = Y, optns = list("lower" = 0)),
               fastfrechet::frechetreg_univar2wass(X, Y, lower = 0),
               times = 5)
```
![**Figure 1**. Fréchet regression optimization accuracy and median solve time
(of 5 iterations) on simulated zero-inflated negative binomial responses, zoomed
in around zero. Departures $\widehat{\mathbf{q}}_i < 0$ are violations of lower
bound.\label{fig:frechetreg_comparison}](figures/frechetreg_comparison.png){width="10in"}

### The Variable Selection Problem

The R package `fastfrechet` solves the variable selection procedure proposed by
@tucker_variable_2023 for Fréchet regression, specifically in 2-Wasserstein
space. The implementation is a second-order geodesic gradient descent algorithm
developed by @coulter_fast_2024, with two new modifications. First, the new
implementation uses the custom dual active-set method discussed in the previous
subsection, with warm starts. Second, the new implementation includes an option
for the user to specify an impulse parameter for the gradient descent algorithm,
which implements a momentum-based geodesic gradient descent.

The algorithm implemented by `fastfrechet` provides a fast and accurate solution
to the variable selection problem, illustrated by comparison to the algorithm
available from @tucker_variable_2023 (Supplementary Material) in Figure
\autoref{fig:friso_comparison}. We use the same $\epsilon = 0.0075$ error
tolerance for the `fastfrechet` method as @coulter_fast_2024, chosen to provide
similar optimization accuracy so computation time comparisons are on
approximately equal footing. Note that `fastfrechet` centers and scales
$\mathbf{X}$ so $\mathrm{diag}(\mathbf{X}^{\top}\mathbf{X}) = \mathbf{1}_p$, to
remove the effect of unit choice on the variable selection outcome. We use this
same centered-scaled $\mathbf{X}$ as input to the old method. Finally, as the
existing code for the old method does not natively accept multiple $\tau$
inputs, we wrote a short wrapper function which loops over the values in
`tauseq`.
```
# Centering and scaling X
X0 = X - rep(1, n) %*% crossprod(rep(1 / n, n), X)
X0 = X0 / (rep(1, n) %*% sqrt(crossprod(rep(1 / n, n), X0 * X0)))

# Defining tau range:
tauseq = seq(0.5, 10, 0.5)

microbenchmark(fastfrechet::FRiSO_univar2wass(X, Y, lower = 0, tauseq = tauseq, eps = 0.0075, nudge = 0.01),
               FRiSO_tucker(X = X0, Y = Y, tauseq = tauseq, lower = 0, upper = 1000),
               times = 5)
```
![**Figure 2**. Variable selection comparison between old method and new method
from `fastfrechet`, across $\tau = \{0.5, 1, \dots, 10\}$. Each algorithm solves
$\widehat{\pmb{\lambda}}(\tau)$ on increasing $\tau$, using warm starts.
(*left*) Variable selection solution paths and computation times. (*right*)
Optimization accuracy comparison; values below zero (grey dotted line)
correspond to superior optimization accuracy of `fastfrechet` method, and values
above zero correspond to superior optimization accuracy of old method.
\label{fig:friso_comparison}](figures/friso_comparison.png){width="8.343585in"}


# Acknowledgements

# References
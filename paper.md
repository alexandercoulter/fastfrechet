---
title: 'fastfrechet: An R package for fast implementation of Fréchet regression with
  distribution responses'
tags:
- Fréchet regression
- Variable selection
- Distributions
- Wasserstein
- Metric space
- R
- C++
date: 7 February 2025
output:
  html_document:
    df_print: paged
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
affiliations:
- name: Department of Statistics, Texas A&M University, United States
  index: 1
- name: Department of Biostatistics, University of Michigan, United States
  index: 2
---

# Summary

Distribution-as-response regression problems are gaining wider attention,
especially within biomedical settings where observation-rich patient specific
data sets are available, such as continuous glucose monitoring
[@matabuena_glucodensities_2021], actigraphy [@ghosal_distributional_2023], and
feature densities in CT scans [@petersen_wasserstein_2021].
To accommodate the complex structure of such problems, @petersen_frechet_2019
proposed a regression framework called *Fréchet regression* which allows
constrained responses. This regression framework was further extended for
variable selection by @tucker_variable_2023, and @coulter_fast_2024 developed a
fast variable selection algorithm for the specific setting of univariate
distribution responses equipped with the 2-Wasserstein metric
(*2-Wasserstein space*). We present `fastfrechet`, an R package providing fast
implementation of these Fréchet regression and variable selection methods in
2-Wasserstein space, with resampling tools for automatic variable selection.
`fastfrechet` makes distribution-based Fréchet regression with
resampling-supplemented variable selection readily available and highly scalable
to large data sets, such as the UK Biobank [@doherty_large_2017].

<!--Distribution-as-response regression problems are gaining wider attention,
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
to large data sets, such as the UK Biobank [@doherty_large_2017].-->

# Statement of Need

No software package currently supports variable selection for 2-Wasserstein
Fréchet regression. Implementation of Fréchet regression in 2-Wasserstein space
without variable selection is supported by two R packages: `WRI` and `frechet`.
Aside from lack of variable selection functionality, these packages face certain
practical limitations. For instance, `WRI` requires continuous distributions,
and does not allow user-specified constraints for the distribution support.
`frechet` offers more flexibility in user specifications, but its solver for
Fréchet regression can be computationally slow and may not satisfy constraints
with sufficient accuracy.

<!--No software package currently supports variable selection for 2-Wasserstein
Fréchet regression. Implementation of Fréchet regression in 2-Wasserstein space without variable selection is supported by two R packages: `WRI` and `frechet`.  Aside from lack of variable selection functionality, these packages face certain practical limitations. For instance, `WRI` requires strictly increasing quantile function inputs and does not allow user-specified constraints for the distribution support. `frechet` offers more flexibility in user specifications but its solver for Fréchet regression can be computationally slow and may not satisfy constraints with sufficient accuracy.-->

The `fastfrechet` package addresses these limitations by providing a fast, scalable, and user-friendly implementation of variable selection for 2-Wasserstein Fréchet regression based on the work of @coulter_fast_2024. The package incorporates resampling tools, including cross-validation as described in @tucker_variable_2023 and stability selection as discussed in @coulter_fast_2024, to enhance automatic variable selection. Additionally, `fastfrechet` features a customized dual active-set solver, inspired by @arnstrom_dual_2022, which ensures both computational efficiency and accuracy while accommodating user-specified constraints in Fréchet regression. The Fréchet regression solver also accommodates the weighting scheme used in the variable selection procedure, the first package to do so.

<!--The `fastfrechet` package addresses these limitations by providing a fast, scalable, and user-friendly implementation of variable selection for 2-Wasserstein Fréchet regression based on the work of @coulter_fast_2024, including the general setting without variable selection as a special case. The package incorporates resampling tools, including cross-validation as described in @tucker_variable_2023 and stability selection as discussed in @coulter_fast_2024, to enhance automatic variable selection. Additionally, `fastfrechet` features a customized dual active-set solver, inspired by @arnstrom_dual_2022, which ensures both computational efficiency and accuracy while accommodating user-specified constraints.-->

<!--The resulting algorithm is  ....[say something about how fast it is so we can apply it now to huge datasets]-->

<!--No software package currently supports variable selection for 2-Wasserstein
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
specification limitations in existing packages.-->

# Performance Comparisons to Existing Implementations

We illustrate the performance of R package `fastfrechet` against existing
implementations with simulated zero-inflated negative binomial distributions.
The included function `generate_zinbinom_qf` allows us to simulate $n$ such
distributions, represented as quantile functions evaluated on a shared $m$-grid
in $(0, 1)$, dependent on the first 4 of $p \geq 4$ covariates. For a general
illustration of what Fréchet regression and variable selection look like in
2-Wasserstein space, see the accompanying `intro-fastfrechet` vignette.
```
library(fastfrechet)

set.seed(31)
n = 100          # number of quantile functions
p = 10           # number of covariates
m = 100          # (0, 1) quantile functions grid density 

gendata = generate_zinbinom_qf(n, p, m)
X = gendata$X    # (n x p) covariate matrix
Y = gendata$Y    # (n x m) quantile response matrix, stored row-wise
```

### The Fréchet Regression Problem

The R package `fastfrechet` provides a solver for the Fréchet regression problem
for 2-Wasserstein space, with optional user-specified  `lower` and `upper` box
constraints to enforce finite support on the underlying distributions. The
implementation is a customization of the dual active-set method of
@arnstrom_dual_2022 (see the accompanying `monotoneQP-fastfrechet` vignette).
The resulting algorithm provides a fast and numerically exact solution to the
Fréchet regression problem, illustrated in \autoref{fig:frechetreg_comparison}.
Since the `WRI` package requires strictly monotone inputs, we add a small
adjustment to the response matrix when using `WRI`'s solver.
```
library(microbenchmark)

X_df = as.data.frame(X)
Y_adj = Y + matrix(seq(0, 0.01, len = m), n, m, byrow = TRUE)

microbenchmark(WRI::wass_regress( ~ ., X_df, "quantile", Y_adj),
               frechet::GloDenReg(X, qin = Y, optns = list("lower" = 0)),
               frechetreg_univar2wass(X, Y, lower = 0),
               times = 5)
```

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
available from @tucker_variable_2023 (Supplementary Material) in 
\autoref{fig:friso_comparison}. We use the same $\epsilon = 0.0075$ error
tolerance for the `fastfrechet` method as @coulter_fast_2024, chosen to provide
similar optimization accuracy so computation time comparisons are on
approximately equal footing. As existing code for the old method does not
natively accept multiple $\tau$ inputs, we wrote a short wrapper function which
loops over the $\tau$ values.
```
# Centering and scaling X
X0 = (scale(X) * sqrt(n / (n - 1)))[,]

# Defining tau range:
tauseq = seq(0.5, 10, 0.5)

microbenchmark(FRiSO_univar2wass(X, Y, lower = 0, tauseq = tauseq,
                                 eps = 0.0075, nudge = 0.01),
               FRiSO_tucker(X0, Y, tauseq = tauseq, lower = 0, upper = 1000),
               times = 5)
```

# Figures

![Fréchet regression optimization accuracy and median solve time
(of 5 iterations) on simulated zero-inflated negative binomial responses, zoomed
in around zero. Departures $\widehat{\mathbf{q}}_i < 0$ are violations of lower
bound.\label{fig:frechetreg_comparison}](figures/frechetreg_comparison.png){width="5.4in"}

![Variable selection comparison between old method and new method
from `fastfrechet`, across $\tau = \{0.5, 1, \dots, 10\}$. Each algorithm solves
$\widehat{\pmb{\lambda}}(\tau)$ on increasing $\tau$, using warm starts.
(*left*) Variable selection solution paths and computation times. (*right*)
Optimization accuracy comparison; values below zero (grey dotted line)
correspond to superior accuracy of `fastfrechet` method, and values
above zero correspond to superior accuracy of old method.
\label{fig:friso_comparison}](figures/friso_comparison.png){width="5.4in"}


# Acknowledgements
IG's research was supported by NSF DMS-2422478 and NIH R01HL172785.

# References
---
title: 'fastfrechet: An R package for fast implementation of Fréchet regression with
  distributional responses'
tags:
- Fréchet regression
- Variable selection
- Distributions
- Wasserstein
- Metric space
- R
- C++
date: 07 March 2025
output:
  html_document:
    df_print: paged
authors:
- name: Alexander Coulter
  orcid: "0009-0001-9116-6588"
  affiliation: 1
  corresponding: yes
- name: Rebecca Lee
  orcid: "0000-0001-8162-6555"
  affiliation: 1
- name: Irina Gaynanova
  orcid: "0000-0002-4116-0268"
  affiliation: 2
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
data sets are available, such as feature densities in CT scans
[@petersen_wasserstein_2021], actigraphy [@ghosal_distributional_2023], and 
continuous glucose monitoring
[@coulter_fast_2024; @matabuena_glucodensities_2021]. To accommodate the
complex structure of such problems, @petersen_frechet_2019 proposed a regression
framework called *Fréchet regression* which allows non-Euclidean responses,
including distributional responses. This regression framework was further
extended for variable selection by @tucker_variable_2023, and @coulter_fast_2024
developed a fast variable selection algorithm for the specific setting of
univariate distributional responses equipped with the 2-Wasserstein metric
(*2-Wasserstein space*). We present `fastfrechet`, an R package providing fast
implementation of these Fréchet regression and variable selection methods in
2-Wasserstein space, with resampling tools for automatic variable selection.
`fastfrechet` makes distribution-based Fréchet regression with
resampling-supplemented variable selection readily available and highly scalable
to large data sets, such as the UK Biobank [@doherty_large_2017].

# Statement of Need

Fréchet regression with variable selection is currently not implemented by any
software package, available only through the Supplementary Material of
@tucker_variable_2023 (hereafter "`Tucker` materials"). As discussed in
@coulter_fast_2024, the implemented algorithm is prohibitively slow.
Implementation of Fréchet regression in 2-Wasserstein space without variable
selection is supported by the `Tucker` materials, and by two R packages: `WRI` [@liu_wri_2022] and `frechet` [@chen_frechet_2023]. These packages face certain
practical limitations. For instance, `WRI` requires continuous distributions,
and does not allow user-specified constraints for the distribution support.
`frechet` offers more flexibility in user specifications, but its solver for
Fréchet regression is slow and may not accurately satisfy constraints.

The `fastfrechet` package addresses these limitations by providing a fast,
scalable, and user-friendly implementation of both Fréchet regression and
variable selection for 2-Wasserstein space, based on the work of
@coulter_fast_2024. The package incorporates resampling tools to enhance
automatic variable selection, including cross-validation described in
@tucker_variable_2023 and stability selection described in @coulter_fast_2024.
Additionally, `fastfrechet` features a customized dual active-set solver for the
Fréchet regression problem, inspired by @arnstrom_dual_2022, which ensures both
computational efficiency and accuracy while accommodating user-specified support
constraints. The Fréchet regression solver also accommodates the auxiliary
weighting scheme used in the variable selection procedure, the first package to
do so.

# Performance Comparisons to Existing Implementations

We illustrate the performance of `fastfrechet` against existing implementations
with simulated covariate-dependent distributional responses. The included
function `generate_zinbinom_qf` simulates $n$ zero-inflated negative binomial
(**zinbinom**) distributions (we choose $n = 100$), represented as quantile
functions evaluated on a shared $m$-grid in $(0, 1)$ (we choose $m = 100$), and
dependent on the first 4 of $p \geq 4$ covariates (we choose $p = 10$). We
utilize the R package `microbenchmark` [@mersmann_microbenchmark_2024] to
calculate run times, and report median times for each method (Fréchet
regression, variable selection) from 15 iterations; all computations were
performed on an Apple M1 Max chip. To replicate the specific simulation and
comparison settings used in this manuscript, see the accompanying `performanceExample-fastfrechet` vignette.

### The Fréchet Regression Problem

`fastfrechet` provides a solver for the Fréchet regression problem for
2-Wasserstein space, with optional `lower` and `upper` support constraints on
the underlying distributions. Since **zinbinom** distributions are non-negative,
we fix `lower = 0` and `upper = Inf` (or some suitably large number, as
applicable). The regression outputs are fitted quantile functions, which should
be monotone non-decreasing and obey support constraints. The `fastfrechet`
implementation is a customization of the dual active-set method of
@arnstrom_dual_2022. (See the accompanying `monotoneQP-fastfrechet` vignette for
full algorithm description.)

\autoref{fig:frechetreg_comparison} illustrates the speed and accuracy of
Fréchet regression implemented in `fastfrechet` against the `WRI`, `frechet`,
and `Tucker` materials implementations. `WRI` does not accept known support
bounds as input, and fitted responses correspondingly violate the zero lower
bound; `frechet` solutions only approximately satisfy the lower bound. The
`Tucker` materials implementation finds numerically accurate solutions, but
`fastfrechet` accomplishes this in a fraction the time.

### The Variable Selection Problem

The R package `fastfrechet` implements variable selection for Fréchet
regression, specifically in 2-Wasserstein space. Variable selection comprises
finding optimal weight vector $\widehat{\pmb{\lambda}} \in \mathbb{R}^p$ that
satisfies a $\tau$-simplex constraint, given hyperparameter $\tau > 0$. In
2-Wasserstein space, $\widehat{\pmb{\lambda}}$ essentially minimizes an $L^2$
norm between weighted Fréchet regression outputs
$\widehat{\pmb{Q}}(\widehat{\pmb{\lambda}})$ and the raw data $\pmb{Y}$. (See
the accompanying `intro-fastfrechet` vignette for a detailed exposition.)
`fastfrechet` implements the second-order geodesic descent algorithm developed
by @coulter_fast_2024, with two modifications. First, the implementation uses
the custom dual active-set method discussed in the previous subsection. The
active set defining the weighted Fréchet regression solution
$\widehat{\pmb{Q}}(\pmb{\lambda}^t)$ for iterate $\pmb{\lambda}^t$ serves as a
warm start for iterate $\pmb{\lambda}^{t + 1}$, reducing computation time.
Second, the implementation allows the user to specify an impulse parameter,
which implements momentum-based geodesic descent.

\autoref{fig:friso_comparison} illustrates the speed and accuracy of variable
selection implemented in `fastfrechet` against the `Tucker` materials
implementation, across sequence of hyperparameter values
$\tau \in \{0.5, 1.0, \cdots, 10.0\}$. We hand-select `fastfrechet` error
tolerance parameter $\varepsilon = 0.014$, which gives solutions
$\widehat{\pmb{\lambda}}(\tau)$ minimizing the objective function
approximately as well as solutions from the other method "as-is". `fastfrechet`
is upward of 20,000$\times$ faster to obtain these comparable solutions.
Decreasing the `fastfrechet` error tolerance parameter increases optimization
accuracy with modest increases in computation time.

# Figures

![Fitted Fréchet regression quantile functions (zoomed in around zero) and
median run times for `fastfrechet` and other implementations. Fitted quantile
functions below zero violate known lower support constraints.\label{fig:frechetreg_comparison}](figures/frechetreg_accuracy_comparison.png){width=100%}

![(**left, center**) Variable selection solution paths
$\widehat{\pmb{\lambda}}(\tau)$ across $\tau = \{0.5, 1, \dots, 10\}$ and median
run times for `Tucker` materials and `fastfrechet`. (**right**) Relative
optimization accuracy of `fastfrechet` and `Tucker` materials variable
selection, and median `fastfrechet` run times, using different error tolerance
values. Points below 1.0 indicate `fastfrechet` solutions minimize the objective
function better.\label{fig:friso_comparison}](figures/friso_accuracy_comparison.png){width=100%}


# Acknowledgements
IG's research was supported by NSF DMS-2422478 and NIH R01HL172785.

# References
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
date: 15 February 2025
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
[@coulter_variable_2024, @matabuena_glucodensities_2021]. To accommodate the
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
to large data sets, such as the UK Biobank [@doherty_large_2017]. For a general
illustration of Fréchet regression and variable selection in 2-Wasserstein
space, see the accompanying `intro-fastfrechet` vignette.

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

Fréchet regression with variable selection is currently not implemented by any
software package, available only through the Supplementary Material of
@tucker_variable_2023 (hereafter `Tucker` materials). As discussed in
@coulter_fast_2024, the implemented algorithm is prohibitively slow.
Implementation of Fréchet regression in 2-Wasserstein space without variable
selection is supported by the `Tucker` materials, and by two R packages: `WRI`
and `frechet`. These packages face certain practical limitations. For instance,
`WRI` requires Lebesgue-continuous distributions, and does not allow
user-specified constraints for the distribution support. `frechet` offers more
flexibility in user specifications, but its solver for Fréchet regression can be
computationally slow and may not satisfy constraints with sufficient accuracy.

<!--No software package currently supports variable selection for 2-Wasserstein
Fréchet regression. Implementation of Fréchet regression in 2-Wasserstein space without variable selection is supported by two R packages: `WRI` and `frechet`.  Aside from lack of variable selection functionality, these packages face certain practical limitations. For instance, `WRI` requires strictly increasing quantile function inputs and does not allow user-specified constraints for the distribution support. `frechet` offers more flexibility in user specifications but its solver for Fréchet regression can be computationally slow and may not satisfy constraints with sufficient accuracy.-->

The `fastfrechet` package addresses these limitations by providing a fast,
scalable, and user-friendly implementation of both Fréchet regression and
variable selection for 2-Wasserstein space, based on the work of
@coulter_fast_2024. The package incorporates resampling tools, including
cross-validation described in @tucker_variable_2023 and stability selection
described in @coulter_fast_2024, to enhance automatic variable selection.
Additionally, `fastfrechet` features a customized dual active-set solver for the
Fréchet regression problem, inspired by @arnstrom_dual_2022, which ensures both
computational efficiency and accuracy while accommodating user-specified support
constraints. The Fréchet regression solver also accommodates the weighting
scheme used in the variable selection procedure, the first package to do so.

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

We illustrate the performance of `fastfrechet` against existing implementations
with simulated covariate-dependent distribution responses. The included function
`generate_zinbinom_qf` allows us to simulate $n$ zero-inflated negative binomial
(**zinbinom**) distributions (we choose $n = 100$), represented as quantile
functions evaluated on a shared $m$-grid in $(0, 1)$ (we choose $m = 100$), and
dependent on the first 4 of $p \geq 4$ covariates (we choose $p = 10$). To
replicate the specific simulation and comparison settings used in this
manuscript, see the accompanying `performanceExample-fastfrechet` vignette.

<!--We illustrate the performance of R package `fastfrechet` against existing
implementations with simulated zero-inflated negative binomial distributions.
The included function `generate_zinbinom_qf` allows us to simulate $n$ such
distributions, represented as quantile functions evaluated on a shared $m$-grid
in $(0, 1)$, dependent on the first 4 of $p \geq 4$ covariates. For a general
illustration of what Fréchet regression and variable selection look like in
2-Wasserstein space, see the accompanying `intro-fastfrechet` vignette.-->

### The Fréchet Regression Problem

`fastfrechet` provides a solver for the Fréchet regression problem for
2-Wasserstein space, with optional user-specified  `lower` and `upper` support
constraints on the underlying distributions. Since **zinbinom** distributions
are non-negative, we fix `lower = 0` and `upper = Inf` (or some suitably large
number). The regression outputs are fitted quantile functions, which should be
monotone non-decreasing and obey support constraints. The `fastfrechet`
implementation is a customization of the dual active-set method of
@arnstrom_dual_2022 (see the accompanying `monotoneQP-fastfrechet` vignette).

\autoref{fig:frechetreg_comparison} illustrates the comparative speed and
accuracy of Fréchet regression implemented in `fastfrechet` against the
implementations from the `WRI` and `frechet` packages, and the `Tucker`
materials. `WRI` does not accept known support bounds as input, and fitted
responses correspondingly violate the zero lower bound; `frechet` solutions only
approximately satisfy the lower bound. Both the `Tucker` and `fastfrechet`
implementations give accurate solutions, but `fastfrechet` obtains the solution
in a fraction the time.

### The Variable Selection Problem

The R package `fastfrechet` implements variable selection for Fréchet
regression, specifically in 2-Wasserstein space. Variable selection comprises
finding an optimal weight vector $\widehat{\pmb{\lambda}} \in \mathbb{R}^p$ that
satisfies a $\tau$-simplex constraint, for fixed hyperparameter $\tau > 0$. In
2-Wasserstein space, $\widehat{\pmb{\lambda}}$ essentially minimizes an $L^2$
norm between weighted Fréchet regression outputs
$\widehat{\pmb{Q}}(\widehat{\pmb{\lambda}})$ and the raw data $\pmb{Y}$. (See
the accompanying `intro-fastfrechet` vignette for an exposition.) `fastfrechet`
implements the second-order geodesic gradient descent algorithm developed by
@coulter_fast_2024, with two new modifications. First, the implementation uses
the custom dual active-set method discussed in the previous subsection, with
warm starts. Second, the implementation includes an option for the user to
specify an impulse parameter, which implements momentum-based geodesic gradient
descent.

\autoref{fig:friso_comparison} illustrates the comparative speed and accuracy of
variable selection implemented in `fastfrechet` against the implementation from
the `Tucker` materials. For a sequence of hyperparameter values
$\tau \in \{0.5, 1.0, \cdots, 10.0\}$, the computation time of `fastfrechet`
is upward of 20,000$\times$ less to obtain solutions
$\widehat{\pmb{\lambda}}(\tau)$ of comparable
objective-minimizing quality, compared to those obtained from implementing the
the `Tucker` materials "as-is". We also show reducing the `fastfrechet` error
tolerance parameter increases optimization accuracy with only modest increases
in computation time.

# Figures

![Fréchet regression optimization accuracy and median solve time
(of 15 iterations) on **zinbinom** responses, zoomed in around zero.\label{fig:frechetreg_comparison}](inst/doc/frechetreg_accuracy_comparison.png){width=100%}

![Variable selection comparison between implementations in `fastfrechet` and the
`Tucker` materials, across $\tau = \{0.5, 1, \dots, 10\}$.\label{fig:friso_comparison}](inst/doc/friso_accuracy_comparison.png){width=100%}


# Acknowledgements
IG's research was supported by NSF DMS-2422478 and NIH R01HL172785.

# References
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

We present `fastfrechet`, an R package providing fast implementation of Fréchet
regression and variable selection methods for univariate distribution responses. 
Distribution-as-response is gaining wider attention, especially within
biomedical settings where observation-rich patient level data sets are
available, such as continuous glucose monitoring
[@matabuena_glucodensities_2021], actigraphy [@ghosal_distributional_2023], and
feature densities in CT scans [@petersen_wasserstein_2021]. Naïve application of
standard Euclidean regression is not appropriate since the response space is not
a vector space. To overcome such issues, @petersen_frechet_2019 proposed an
extension of regression with Euclidean covariates, for responses within a
general metric space, so-called *Fréchet regression*. Using this general
framework, @tucker_variable_2023 proposed a variable selection method, extending
@wu_cant_2021.

The solution to the Fréchet regression problem is specific to the metric space.
As such, the original variable selection implementation employs a coordinate
descent algorithm which favors generality with the cost of prohibitively slow
computation time. @coulter_fast_2024 developed a fast algorithm to solve the
variable selection problem for the specific setting of univariate distributions
equipped with the 2-Wasserstein metric (*2-Wasserstein space*). `fastfrechet`
implements this new variable selection solver, along with a customized dual
active-set method [@arnstrom_dual_2022] for solving the Fréchet regression
problem efficiently. It also includes implementations of the cross-validation
and stability selection procedures respectively described in
@tucker_variable_2023 and @coulter_fast_2024. These improvements allow Fréchet
regression with resampling-supplemented variable selection to be viable and
readily available procedures for distribution responses, especially for
application to large data sets like the UK Biobank [@doherty_large_2017].

# Statement of Need

For Fréchet regression, two publicly available software packages exist which
implement Fréchet regression in 2-Wasserstein space. The R package `WRI` has
several practical limitations: it requires strictly monotone empirical quantile
functions (or, strictly positive empirical densities), which is unnecessary in
this setting; it sets a seemingly arbitrary minimum threshold of 25 observations
per empirical response; and the user cannot specify box constraints for bounded
distribution support, so fitted responses can violate known bounds even when
observed responses do not. By contrast, the R package `frechet` largely does not
have these limitations, however its solver still does not reliably obey box
constraints and is relatively slow. The R package `fastfrechet` does not have
the practical limitations of `WRI`, while it obtains more accurate solutions in
a fraction the time compared to `frechet`.

For variable selection, the coordinate descent algorithm is available as
part of the supplementary material to @tucker_variable_2023. However, this
implementation is not readily accessible through a repository, nor in package
structure. Part of the variable selection problem is solving a weighted version
of the associated Fréchet regression problem, but neither `WRI` nor `frechet` is
amenable to specifying these weights. To our knowledge, there is no existing
package which implements the variable selection procedure of
@tucker_variable_2023, either fully or partially. The R package `fastfrechet`
provides ready implementation of the variable selection method through a fast,
robust gradient descent algorithm, and allows the user to optionally specify the
weighting structure for Fréchet regression in the variable selection framework.

# Performance Comparisons

We illustrate the performance of R package `fastfrechet` against existing
implementations with simulated covariate-dependent, zero-inflated negative
binomial distributions. The included function `generate_zinbinom_qf` allows us
to simulate $n$ such distributions, represented as quantile functions evaluated
on a shared $m$-grid in $(0, 1)$, generated from the first 4 of $p \geq 4$
covariates. For a general illustration of what Fréchet regression and variable
selection look like in 2-Wasserstein space, see the accompanying `VIGNETTE_NAME`
vignette.
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
accompanying `VIGNETTE_NAME` vignette). The resulting algorithm provides a fast
and numerically exact solution to the Fréchet regression problem, illustrated in
Figure \autoref{fig:frechetreg_comparison}. Note the `WRI` package requires
strictly monotone inputs, so add a small monotone adjustment to the response
matrix when using that package.
```
library(microbenchmark)

X_df = as.data.frame(X)
Y_adj = Y + matrix(seq(0, 0.01, len = m), n, m, byrow = TRUE)

microbenchmark(WRI::wass_regress( ~ ., X_df, "quantile", Y_adj),
               frechet::GloDenReg(xin = X, qin = Y, optns = list("lower" = 0)),
               fastfrechet::frechetreg_univar2wass(X, Y, lower = 0),
               times = 5)
```
![**Figure 1**. Fréchet regression optimization accuracy and median solve time
(of 5 iterations) on simulated zero-inflated negative binomial responses, zoomed
in around zero. Departures $\widehat{\mathbf{q}}_i < 0$ are violations of lower
bound.\label{fig:frechetreg_comparison}](figures/frechetreg_comparison.png)

### The Variable Selection Problem

The 

# Old writing

Fréchet regression (@petersen_frechet_2019) extends Euclidean regression to the
general setting where response $Y$ resides in a metric space $\Omega$ equipped
with a metric $d : \Omega \times \Omega \mapsto \mathbb{R}_+$. Subsequently, variable
selection has been proposed for Fréchet regression (@tucker_variable_2023),
again extending a generalized ridge penalty method developed in the Euclidean
setting (@wu_cant_2021). Finding conditional Fréchet means (i.e. solving the
Fréchet regression problem) is a context-specific task where general solvers are
not available; on the other hand, @tucker_variable_2023 propose a general method
for solving their variable selection problem, so long as some method for solving
the associated Fréchet regression problem is available.

Univariate distribution responses, equipped with the 2-Wasserstein metric, serve
as an example use case in both @petersen_frechet_2019 and @tucker_variable_2023,
and have emerging applications to biomedical data including continuous glucose
monitoring data [MATABUENA REFERENCE] and actigraphy data [GHOSAL REFERENCE].
@coulter_fast_2024 developed a new algorithm for solving the variable selection
problem in this context, obtaining empirical speed increases up to
10,000$\times$ in simulation settings over the existing general purpose solver,
with greater relative speed increases with increasing sample size and covariate
count. This makes Fréchet regression methods accessible for data sets where such
data is available for a large number of patient and covariates, such as the UK
Biobank [DOHERTY REFERENCE].

`fastfrechet` implements the fast variable selection algorithm of
@coulter_fast_2024 for univariate distribution responses, as well as implements
resampling procedures for variable selection like cross validation
(@tucker_variable_2023) and complementary pairs stability selection
(@coulter_fast_2024; @shah_variable_2013; @meinshausen_stability_2010). It also
includes a new dedicated QP solver for the associated Fréchet regression
problem, implementing the dual active-set method of @arnstrom_dual_2022 while
taking advantage of the specific constraint structure to avoid matrix
decomposition and multiplication operations.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

Example citation [@coulter_fast_2024].

# Figures

# Acknowledgements

# References
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
date: "12 November 2024"
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
available, such as continuous glucose monitoring [MATABUENA REFERNECE],
actigraphy [GHOSAL REFERENCE], and feature densities in CT scans [PETERSEN, LIU,
DIVANI REFERENCE]. Naïve application of standard Euclidean regression is not
appropriate since the response space is not a vector space. To overcome such
issues, @petersen_frechet_2019 proposed an extension of regression with
Euclidean covariates, for responses within a general metric space, so-called
*Fréchet regression*. Using this general framework, @tucker_variable_2023
proposed a variable selection method through a similar extension of a Euclidean
method (@wu_cant_2021). As the solution to the Fréchet regression problem is
specific to the metric space, the initially published implementation of the
variable selection procedure employs a coordinate descent algorithm which makes
minimal use of the problem geometry, favoring generality at the expense of
prohibitively slow computation time. @coulter_fast_2024 developed a second order
gradient descent algorithm for the specific setting of univariate distributions
equipped with the 2-Wasserstein metric (*2-Wasserstein space*), decreasing
solution time on the order of 10,000\eqn{\times}. `fastfrechet` implements this
new variable selection solver, along with a customized active set method
(@arnstrom_dual_2022) for solving the Fréchet regression problem efficiently. It
also includes implementations of the cross-validation and stability selection
procedures respectively described in @tucker_variable_2023 and
@coulter_fast_2024, allowing Fréchet regression with resampling-supplemented
variable selection to be a viable and readily available procedure for
distribution response settings.

# Statement of Need

For Fréchet regression, two publicly available software package exists which
implement Fréchet regression in 2-Wasserstein space. The R package `WRI` has
several practical limitations: it requires strictly monotone empirical quantile
functions (or, strictly positive empirical densities), which is unnecessary in
this setting; it sets a seemingly arbitrarily minimum threshold of 25
observations per empirical response; and the user cannot specify box constraints
for finite distribution support, so fitted responses can violate known bounds
even when observed responses do not. By contrast, the R package `frechet`
largely does not have these limitations, however its solver still does not
reliably obey box constraints and is relatively slow. The R package
`fastfrechet` does not have the practical limitations of `WRI`, while it obtains
more accurate solutions in a fraction the time compared to `frechet`.

For variable selection, the coordinate descent algorithm is available as
part of the supplementary material to @tucker_variable_2023. However, this
implementation is not readily accessible through a repository nor in package
structure. Part of the variable selection problem is solving a weighted version
of the associated Fréchet regression problem, but neither `WRI` nor `frechet` is
amenable to specifying these weights. As such, to our knowledge there is no
existing package which implements the variable selection procedure of
@tucker_variable_2023, either fully or partially. The R package `fastfrechet`
provides ready implementation of the variable selection method through a fast,
robust gradient descent algorithm, and allows the user to optionally specify the
weighting structure for Fréchet regression in the variable selection framework.

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
problem in this context, obtaining empirical speed increases up to 10,000x in
simulation settings over the existing general purpose solver, with greater
relative speed increases with increasing sample size and covariate count. This
makes Fréchet regression methods accessible for data sets where such data is
available for a large number of patient and covariates, such as the UK Biobank
[DOHERTY REFERENCE].

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
---
title: 'fastfrechet: An R package for fast implementation of Fréchet regression functions
  for distribution responses'
tags:
- R
- Fréchet regression
- Variable selection
- Distributions
- Wasserstein
- Metric space
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

# Statement of need

This is the paper's statement of need [edit from above?].

# Mathematics

This is an equation example $x = 1$.

This is another equation example
$$
\begin{align}
x &= 1 \\
x + 2 &= 3
\end{align}
$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastfrechet: Fast Fréchet regression with univariate distributional responses

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexandercoulter/fastfrechet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexandercoulter/fastfrechet/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/alexandercoulter/fastfrechet/graph/badge.svg)](https://app.codecov.io/gh/alexandercoulter/fastfrechet)
<!-- badges: end -->

`fastfrechet` is an `R` package for performing Fréchet regression with
univariate distributional responses, including variable selection. The
corresponding methodological frameworks have been described in

- [Petersen A and Müller H-G (2019), “Fréchet regression for random
  objects with Euclidean predictors”, *Annals of
  Statistics*](https://doi.org/10.1214/17-AOS1624)
- [Tucker DC, Wu Y and Müller H-G (2023), “Variable Selection for Global
  Fréchet Regression”,
  *JASA*](https://doi.org/10.1080/01621459.2021.1969240)

The fast implementation and stability selection is based on the
algorithm developed in

- [Coulter A, Aurora RN, Punjabi N and Gaynanova I (2024), “Fast
  variable selection for distributional regression with application to
  continuous glucose monitoring data.”
  *arXiv*](https://arxiv.org/abs/2403.00922)

## Installation

To use `fastfrechet`, you need to install
[`R`](https://cran.r-project.org/). To enhance your user experience, you
may use some IDE for it (e.g. [`RStudio`](https://www.rstudio.com/)).

The development version of
[`fastfrechet`](https://github.com/alexandercoulter/fastfrechet) is
available on GitHub. You can download it with the help of the `devtools`
package in `R` as follow:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/alexandercoulter/fastfrechet", build_vignettes = TRUE)
```

## Community Guidelines

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution
    guidelines](.github/CONTRIBUTING.md) prior to submitting a pull
    request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/alexandercoulter/fastfrechet/issues).
3.  Contributors must adhere to the [Code of
    Conduct](CODE_OF_CONDUCT.md).

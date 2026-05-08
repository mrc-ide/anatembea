# anatembea

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/anatembea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/anatembea/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

*anatembea* is a tool for AnteNatal Assessment of Temporal malaria Epidemiology using Mechanistic models, Bayesian Estimation, and Analysis, an open-source R package, which translates monthly malaria prevalence collected at first antenatal clinic visit (ANC1) into seasonal patterns in malaria transmission and clinical cases. Our framework is underpinned by an established mechanistic model of malaria dynamics, providing the basis for simulating data-informed hypothetical scenarios to evaluate the likely impact of alternative intervention approaches.

For more information and detail please see our pre-print: <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5184997>

## Installation

You can install the development version of *anatembea* from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/anatembea")
library(anatembea)
```

## Reproducibility for accepted ANC pMCMC manuscript

This package version has been frozen to preserve the exact state used in the accepted manuscript, published in **The Lancet Microbe**. The accepted manuscript is available from the Spiral Digital Repository: <https://hdl.handle.net/10044/1/128591>

The analyses were run with R 4.4 using this frozen package state. To install the paper version, use:

``` r
remotes::install_github("mrc-ide/anatembea@v1.0-paper")
```

Key package versions used in the manuscript workflow included:

- `odin` 1.5.11
- `dust` 0.15.3
- `mcstate` 0.9.22

Analysis scripts and workflow are available at:
<https://github.com/jt-hicks/anc_pmcmc_pub>

(Full citation will be added once the article is published online.)

## Usage

For an example implementation, check out the at <https://mrc-ide.github.io/anatembea/articles/anatembea_intro.html>

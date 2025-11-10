# Compare function to calculate likelihood

`compare_mg` Compare function that compares observed data with model
estimate to calculate likelihood for the particle filter. Converts under
5 year old prevalence in the model to prevalence among multigravid women
based on coefficients from a separate regression analysis.

## Usage

``` r
compare_mg(state, observed, pars = NULL)
```

## Arguments

- state:

  Model output. Default = NULL

- observed:

  Oberved data. Default = NULL

- pars:

  Parameters, optional.

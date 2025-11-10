# Compare function to calculate likelihood

`compare_u5` Compare function that compares observed data with model
estimate to calculate likelihood for the particle filter. Equates
observed prevalence to prevalence under 5 year olds in the model.

## Usage

``` r
compare_u5(state, observed, pars = NULL)
```

## Arguments

- state:

  Model output. Default = NULL

- observed:

  Oberved data. Default = NULL

- pars:

  Parameters, optional.

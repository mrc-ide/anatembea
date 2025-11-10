# Compare function to calculate likelihood

`compare_ancall` Compare function that compares observed data with model
estimate to calculate likelihood for the particle filter. Converts under
5 year old prevalence in the model to prevalence among all women based
on coefficients from a separate regression analysis (TZ data only).

## Usage

``` r
compare_ancall(state, observed, pars = NULL)
```

## Arguments

- state:

  Model output. Default = NULL

- observed:

  Oberved data. Default = NULL

- pars:

  Parameters, optional.

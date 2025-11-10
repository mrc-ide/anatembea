# Compare function to calculate likelihood when fitting to a flexible age group

`compare_flex` Compare function that compares observed data with model
estimate to calculate likelihood for the particle filter. For fitting to
a general population of a user-specified age group

## Usage

``` r
compare_flex(state, observed, pars = NULL)
```

## Arguments

- state:

  Model output. Default = NULL

- observed:

  Oberved data. Default = NULL

- pars:

  Parameters, optional.

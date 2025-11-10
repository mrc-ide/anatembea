# Function that returns values from optional seasonal deterministic model

`check_seasonality` Calculates time series of the prefix seasonal
deterministic model given the posterior distribution of pMCMC parameters

## Usage

``` r
check_seasonality(theta, mpl_pf, season_model)
```

## Arguments

- theta:

  Posterior distribution of pMCMC parameters. Default = NULL

- mpl_pf:

  Model parameter list. Default = NULL

- season_model:

  Seasonality model to be used for the optional deterministic model.
  Default = NULL

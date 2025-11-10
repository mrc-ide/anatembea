# Transformation function that calculates initial values for stochastic model

`data_informed` Calculates the model equilibrium based on an initial EIR
values, then optionally runs a deterministic seasonal model and returns
initial values to be used for the stochastic model fitting.

## Usage

``` r
data_informed(mpl_pf, season_model)
```

## Arguments

- mpl_pf:

  Model parameter list. Default = NULL

- season_model:

  Seasonality model to be used for the optional deterministic model.
  Default = NULL

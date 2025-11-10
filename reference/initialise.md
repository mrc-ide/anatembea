# Estimate the initial state given user inputs

`initialise` Calculates the model equilibrium based on an initial EIR
values, then optionally runs a deterministic seasonal model and returns
initial values to be used for the stochastic model fitting.

## Usage

``` r
initialise(init_EIR, mpl, det_model)
```

## Arguments

- init_EIR:

  User supplied initial EIR

- mpl:

  Model parameter list. Default = NULL

- det_model:

  Seasonality model to be used for the optional deterministic model.
  Default = NULL

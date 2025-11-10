# Equilibrium initialisation list creation

`equilibrium_init_create_stripped` creates an equilibrium initialisation
state to be used within later model runs

## Usage

``` r
equilibrium_init_create_stripped(
  age_vector,
  het_brackets,
  country = NULL,
  admin_unit = NULL,
  ft,
  init_EIR,
  model_param_list,
  state_check = 0
)
```

## Arguments

- age_vector:

  Vector of age brackets.

- het_brackets:

  Integer number of biting heteogenity compartments.

- country:

  String for country of interest. If NULL the seasonal parameters will
  attempt to be loaded using just the admin unit, however if there is
  ambiguity in the admin unit an error will be thrown. If both NULL then
  no seasonality is assumed. Default = NULL.

- admin_unit:

  String for admin unit with country for loading seasonal parameters. If
  country is NULL, the admin unit will attempt to be located,however if
  there is ambiguity in the admin unit an error will be thrown. If both
  country and admin_unit are NULL then no seasonality is assumed.
  Default = NULL.

- ft:

  Numeric for the frequency of people seeking treatment.

- init_EIR:

  Numeric for desired annual EIR.

- model_param_list:

  List of epidemiological parameters created by

- state_check:

  If state_check = TRUE, returns expected deriv values which should
  equal 0 and sets stochastic model to have EIR constant at init_EIR If
  state_check = TRUE and seasonality_on = 1, then the deterministic
  seasonal model is still run, but theta2 is forced to 1, forcing a
  constant seasonality profile If state_check = FALSE, no values are
  printed

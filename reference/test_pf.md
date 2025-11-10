# Test the particle filter

This function sets up and runs a particle MCMC that uses Dust, Odin and
MCState

## Usage

``` r
test_pf(
  data_raw = NULL,
  data_raw_pg = NULL,
  data_raw_mg = NULL,
  n_particles = 200,
  proposal_matrix,
  max_param = 1000,
  max_steps = 1e+07,
  atol = 1e-06,
  rtol = 1e-06,
  n_steps = 500,
  n_threads = 4,
  lag_rates = 10,
  state_check = 0,
  country = NULL,
  admin_unit = NULL,
  preyears = 2,
  seasonality_on = 1,
  seasonality_check = 0,
  seed = 1L,
  start_pf_time = 30,
  stoch_param = c("EIR", "betaa"),
  comparison = c("u5", "pgmg")
)
```

## Arguments

- data_raw:

  Time series data to fit model

- data_raw_pg:

  Time series of primigravidae ANC data to fit model

- data_raw_mg:

  Time series of multigravidae ANC data to fit model

- n_particles:

  Number of particles to be used in pMCMC (default = 200)

- proposal_matrix:

  Proposal matrix for MCMC parameters

- max_param:

  Ceiling for proposed stochastic parameter (either EIR or betaa) values
  (default = 1000)

- max_steps:

  Maximum steps for particle filter (default = 1e7)

- atol:

  atol for particle filter (default = 1e-3)

- rtol:

  rtol for particle filter (default = 1e-6)

- n_steps:

  Number of MCMC steps in a single chain (default = 500)

- n_threads:

  Number of processing threads (default = 4)

- lag_rates:

  Number of delay compartments (default = 10)

- state_check:

  Run equilibrium checks, if state_check = 1, returns expected deriv
  values which should equal 0 and sets stochastic model to have EIR
  constant at init_EIR If state_check = 1 and seasonality_on = 1, then
  the deterministic seasonal model is still run, but theta2 is forced to
  1, forcing a constant seasonality profile If state_check = 0, no
  values are printed

- country:

  Name of country (needed if using seasonality model)

- admin_unit:

  Name of administrative unit (needed if using seasonality model)

- preyears:

  Length of time in years the deterministic seasonal model should run
  before Jan 1 of the year observations began (default = 2)

- seasonality_on:

  Toggle seasonality model run before observed time period (default = 1)

- seasonality_check:

  Toggle saving values of seasonality equilibrium (default = 1)

- seed:

  Allows user to specify a seed (default = 1L)

- start_pf_time:

  Number of days before first observation that particle filter will
  start (default = 30)

- stoch_param:

  Which parameter to test

- comparison:

  The comparison function to be used. Either 'u5' which equates the
  observed prevalence to prevalence under 5 years old in the model or
  'pgmg' which calculates prevalence in primigravid and multigravid
  pregnant women for comparison with observed ANC data.
  c('u5','pg','sg','mg','pgmg','pgsg','ancall')

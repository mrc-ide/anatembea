# Run a pMCMC

This function sets up and runs a particle MCMC that uses Dust, Odin and
MCState

## Usage

``` r
run_pmcmc(
  data_raw = NULL,
  data_raw_pg = NULL,
  data_raw_mg = NULL,
  init_EIR = NULL,
  target_prev = NULL,
  target_prev_group = "u5",
  n_particles = 200,
  proposal_matrix = matrix(1),
  max_param = 125,
  prop_treated = 0.4,
  n_steps = 500,
  n_threads = 4,
  n_chains = 1,
  n_workers = 1,
  state_check = 0,
  country = NULL,
  admin_unit = NULL,
  seasonality_on = TRUE,
  seasonality_check = FALSE,
  check_flexibility = FALSE,
  seed = 1L,
  start_pf_time = 30 * 12,
  particle_tune = FALSE,
  comparison = "u5",
  initial = "informed"
)
```

## Arguments

- data_raw:

  Time series data to fit model

- data_raw_pg:

  Time series of primigravidae ANC data to fit model

- data_raw_mg:

  Time series of multigravidae ANC data to fit model

- init_EIR:

  A single value or a dataframe with two columns (time and EIR) to
  specify historical malaria transmission levels before data collection
  began.

- target_prev:

  Return an initial EIR value (from the equilibrium solution), given a
  target prevalence in under 5yos or 2 to 10 year olds

- target_prev_group:

  Age group used for target prevalence ('u5' or '2to10')

- n_particles:

  Number of particles to be used in pMCMC (default = 200)

- proposal_matrix:

  Proposal matrix for MCMC parameters

- max_param:

  Ceiling for proposed stochastic parameter (either EIR or betaa) values
  (default = 1000)

- prop_treated:

  Proportion of clinical cases that receive effective treatment (default
  = 40%)

- n_steps:

  Number of MCMC steps in a single chain (default = 500)

- n_threads:

  Number of processing threads (default = 4)

- n_chains:

  Number of chains (default = 1)

- n_workers:

  Number of workers (default = 4)

- state_check:

  If state_check = TRUE, returns expected deriv values which should
  equal 0 and sets stochastic model to have EIR constant at init_EIR If
  state_check = TRUE and seasonality_on = 1, then the deterministic
  seasonal model is still run, but theta2 is forced to 1, forcing a
  constant seasonality profile If state_check = FALSE, no values are
  printed

- country:

  Name of country (needed if using seasonality model)

- admin_unit:

  Name of administrative unit (needed if using seasonality model)

- seasonality_on:

  Toggle seasonality model run before observed time period (default = 1)

- seasonality_check:

  Toggle saving values of seasonality equilibrium (default = 1)

- check_flexibility:

  Toggle saving values of flexibility period

- seed:

  Allows user to specify a seed (default = 1L)

- start_pf_time:

  Number of days before first observation that particle filter will
  start (default = 30)

- particle_tune:

  Logical to determine if tuning the number of particles should be
  performed.

- comparison:

  The comparison function to be used. Either 'u5' which equates the
  observed prevalence to prevalence under 5 years old in the model or
  'pgmg' which calculates prevalence in primigravid and multigravid
  pregnant women for comparison with observed ANC data.
  c('u5','pg','sg','mg','pgmg','pgsg','ancall') If in a format XtoY,
  where X and Y are two numbers, will compare to general population
  between those two ages.

- initial:

  Is the initial equilibrium state informed by the user
  ('user-informed') or by the observed data ('fitted')?

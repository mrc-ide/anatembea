# Tuning an anatembea pMCMC

Tuning answers two practical questions: how many particles are needed
for a stable likelihood estimate, and how should the proposal move
through the parameter space? It is a separate stage from scientific
inference.

This article assumes that the data and scientific settings have passed
the checks in [Preparing data and configuration for
anatembea](https://mrc-ide.github.io/anatembea/articles/preparing-an-analysis.md).
It ends with a provisional particle count and proposal covariance. Those
settings are tested with independent chains in [Running and assessing
production pMCMC
chains](https://mrc-ide.github.io/anatembea/articles/production-pmcmc.md).

The numerical values below are starting points, not universal
thresholds. Data length, prevalence, missingness, model options, and
computing resources all affect pMCMC performance.

## The tuning sequence

| Stage | Typical starting design | Question answered |
|----|----|----|
| Smoke test | 100 particles, 10–20 iterations | Does the complete execution path work? |
| Particle check | Repeated filters at 50–400 particles | Is the likelihood stable enough? |
| Covariance pilot | 500 iterations, one chain | Which directions should the proposal move? |
| Scale sweep | Three 500-iteration trials | How far should it move? |

Do not use smoke-test or tuning draws for final inference.

## 1. Build the starting proposal

A proposal covariance describes both the size and direction of proposed
MCMC moves. For informed initialization, only `volatility` is fitted, so
the identity covariance is one-dimensional:

``` r

proposal_informed <- diag(1)
dimnames(proposal_informed) <- list("volatility", "volatility")
proposal_informed
#>            volatility
#> volatility          1
```

For fitted initialization, the order is `log_init_EIR`, then
`volatility`:

``` r

parameter_names <- c("log_init_EIR", "volatility")
proposal_fitted <- diag(length(parameter_names))
dimnames(proposal_fitted) <- list(parameter_names, parameter_names)
proposal_fitted
#>              log_init_EIR volatility
#> log_init_EIR            1          0
#> volatility              0          1
```

The diagonal entries are proposal variances. Off-diagonal entries allow
the parameters to move together. An identity matrix is a neutral way to
test the code path, but it usually ignores the posterior’s scale and
correlation and is therefore not a production proposal.

Always check the matrix before starting a costly run:

``` r

check_covariance <- function(x, parameter_names) {
  stopifnot(
    is.matrix(x),
    identical(dim(x), rep(length(parameter_names), 2L)),
    isTRUE(all.equal(x, t(x))),
    all(is.finite(x))
  )
  if (min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("The proposal covariance must be positive definite.")
  }
  invisible(x)
}

check_covariance(proposal_fitted, parameter_names)
```

## 2. Run a smoke test

The smoke test verifies compilation, arguments, output structure,
versions, and cluster wiring. The example is not evaluated while
building the vignette because model compilation and pMCMC are
deliberately expensive.

``` r

smoke <- anatembea::run_pmcmc(
  data_raw = analysis_data,
  comparison = "u5",
  initial = "fitted",
  check_flexibility = TRUE,
  proposal_matrix = proposal_fitted,
  n_particles = 100,
  n_steps = 20,
  n_threads = 1,
  n_chains = 1,
  n_workers = 1,
  seed = 1101L,
  output_level = "minimal",
  save_state = FALSE,
  save_trajectories = FALSE
)

names(smoke)
head(smoke$pars)
```

`output_level = "minimal"` retains the likelihood inputs but avoids
retaining unneeded model outputs. With `save_trajectories = FALSE`,
`history` and `times` will be `NULL`. This is appropriate for short
tuning runs that only need parameter and likelihood diagnostics.

A successful smoke test is evidence that the workflow executes—not that
the chain converged or that the proposal is suitable.

## 3. Check the particle count

The particle filter produces a noisy estimate of the likelihood. Too few
particles can cause highly variable likelihood estimates, repeated MCMC
states, or filter failures; too many increase runtime without
necessarily improving efficiency.

At a representative fixed parameter vector, repeat the particle filter
with independent seeds at particle counts such as 50, 100, 200, and 400.
Record:

- the variance or standard deviation of the estimated log likelihood;
- particle-filter and ODE failures;
- elapsed time and memory; and
- the exact parameter vector and seed.

The current high-level
[`run_pmcmc()`](https://mrc-ide.github.io/anatembea/reference/run_pmcmc.md)
interface performs a complete pMCMC rather than exposing a standalone
fixed-parameter likelihood replicate. Consequently, this check requires
a project-specific particle-filter wrapper built from the same model
specification. Do not treat changing likelihood values from a single
short MCMC chain as independent likelihood replicates.

A useful results layout is:

``` r

particle_check <- data.frame(
  particles = c(50, 100, 200, 400),
  likelihood_sd = NA_real_,
  failures = NA_integer_,
  elapsed_seconds = NA_real_,
  memory_mb = NA_real_
)

particle_check
#>   particles likelihood_sd failures elapsed_seconds memory_mb
#> 1        50            NA       NA              NA        NA
#> 2       100            NA       NA              NA        NA
#> 3       200            NA       NA              NA        NA
#> 4       400            NA       NA              NA        NA
```

Select the smallest count that provides a sufficiently stable likelihood
with no systematic failures. One hundred particles is often a useful
tuning start and 200 a confirmation or production start, but the
empirical check takes precedence. Hold particle count fixed while
comparing proposal scales.

## 4. Estimate a covariance from a pilot

Run a separate pilot for every materially different time series or
analysis unit. A typical starting design is one 500-iteration chain
using the selected tuning particle count and a unique seed.

``` r

pilot <- anatembea::run_pmcmc(
  data_raw = analysis_data,
  comparison = "u5",
  initial = "fitted",
  check_flexibility = TRUE,
  proposal_matrix = proposal_fitted,
  n_particles = 100,
  n_steps = 500,
  seed = 2101L,
  output_level = "minimal",
  save_state = FALSE,
  save_trajectories = FALSE
)
```

Discard burn-in, retain finite parameter draws, and estimate their
covariance:

``` r

burn_in <- floor(nrow(pilot$pars) / 2)
retained <- pilot$pars[(burn_in + 1):nrow(pilot$pars), parameter_names,
                       drop = FALSE]
retained <- retained[stats::complete.cases(retained), , drop = FALSE]

if (nrow(retained) < 50) {
  stop("Too few usable pilot draws; diagnose or extend the pilot.")
}

pilot_covariance <- stats::cov(retained)
pilot_covariance <- (pilot_covariance + t(pilot_covariance)) / 2

# Correct only small numerical defects; do not hide an unidentified parameter.
eigen_floor <- 1e-8
decomp <- eigen(pilot_covariance, symmetric = TRUE)
if (any(decomp$values <= 0)) {
  stop(
    "The pilot covariance has a non-positive direction; obtain more moving draws."
  )
}
decomp$values <- pmax(decomp$values, eigen_floor)
pilot_covariance <- decomp$vectors %*%
  diag(decomp$values) %*% t(decomp$vectors)
dimnames(pilot_covariance) <- list(parameter_names, parameter_names)

check_covariance(pilot_covariance, parameter_names)
pilot_covariance
```

Review trace plots, parameter ranges, acceptance, repeated states,
failures, and correlation before using this matrix. If the pilot barely
moves or yields fewer than about 50 usable retained draws, do not
conceal the problem with aggressive regularization. Diagnose it or
repeat a longer pilot.

## 5. Sweep the proposal scale

Use the pilot acceptance rate to choose an initial scale grid:

| Pilot acceptance          | Covariance multipliers |
|---------------------------|------------------------|
| Below 0.35 or unavailable | 0.75, 1.0, 1.5         |
| 0.35–0.45                 | 1.0, 1.5, 2.0          |
| Above 0.45                | 1.5, 2.0, 2.83         |

The multiplier here is applied directly to the covariance. If you
instead intend to multiply proposal standard deviations, square that
multiplier before applying it to the covariance.

``` r

pilot_mcmc <- coda::as.mcmc(pilot$pars[, parameter_names, drop = FALSE])
pilot_acceptance <- 1 - coda::rejectionRate(pilot_mcmc)
pilot_acceptance <- min(pilot_acceptance)

scales <- if (pilot_acceptance < 0.35) {
  c(0.75, 1, 1.5)
} else if (pilot_acceptance <= 0.45) {
  c(1, 1.5, 2)
} else {
  c(1.5, 2, 2.83)
}

scale_runs <- Map(function(multiplier, seed) {
  anatembea::run_pmcmc(
    data_raw = analysis_data,
    comparison = "u5",
    initial = "fitted",
    check_flexibility = TRUE,
    proposal_matrix = multiplier * pilot_covariance,
    n_particles = 100,
    n_steps = 500,
    seed = seed,
    output_level = "minimal",
    save_state = FALSE,
    save_trajectories = FALSE
  )
}, scales, c(3101L, 3102L, 3103L))
```

Summarize every parameter, not just the best-looking one:

``` r

summarize_run <- function(result, scale) {
  chain <- coda::as.mcmc(result$pars[, parameter_names, drop = FALSE])
  rejection_rate <- coda::rejectionRate(chain)
  acceptance <- 1 - rejection_rate
  ess <- coda::effectiveSize(chain)
  hours <- as.numeric(result$run_time, units = "hours")

  data.frame(
    scale = scale,
    parameter = names(ess),
    acceptance = acceptance[names(ess)],
    rejection_rate = rejection_rate[names(ess)],
    ess = as.numeric(ess),
    ess_per_hour = as.numeric(ess) / hours
  )
}

scale_summary <- do.call(
  rbind,
  Map(summarize_run, scale_runs, scales)
)
scale_summary
```

Rank successful trials primarily by their minimum ESS per hour across
fitted parameters. Acceptance around 0.25–0.45 is a useful guide for
this low-dimensional random-walk proposal, but trace movement and
absence of failures matter more than hitting a particular number. R-hat
is not informative for a single short chain.

The minimal scale trials do not contain fitted trajectories. After
selecting a provisional multiplier from their parameter diagnostics, run
one explicit trajectory check before confirmation. Here `selected_scale`
is the multiplier chosen from `scale_summary`.

``` r

provisional_check <- anatembea::run_pmcmc(
  data_raw = analysis_data,
  comparison = "u5",
  initial = "fitted",
  check_flexibility = TRUE,
  proposal_matrix = selected_scale * pilot_covariance,
  n_particles = 100,
  n_steps = 500,
  seed = 3201L,
  output_level = "standard",
  save_state = FALSE,
  save_trajectories = TRUE
)

stopifnot(
  !is.null(provisional_check$history),
  !is.null(provisional_check$times)
)
```

Review the fitted prevalence, EIR, emergence, and incidence trajectories
for finite values and scientific plausibility. If every scale performs
poorly, make one targeted adjustment: improve the covariance estimate,
alter the grid, increase particles if the likelihood is noisy, or
investigate model stiffness. Repeatedly increasing chain length without
a diagnosis is not tuning.

## Gate: ready for confirmation?

Freeze a provisional configuration only when you have recorded:

- the particle-count evidence;
- the finite, symmetric, positive-definite covariance and parameter
  order;
- the selected covariance multiplier;
- acceptance, rejection rates, ESS/hour, traces, and failures for every
  trial;
- the provisional fitted-trajectory check;
- all seeds and runtimes; and
- a reason for selecting the winning configuration.

Next, test that configuration with [independent confirmation
chains](https://mrc-ide.github.io/anatembea/articles/production-pmcmc.md).

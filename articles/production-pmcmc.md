# Running and assessing production pMCMC chains

Production begins only after the analysis specification is frozen and
tuning has identified a defensible particle count and proposal
covariance. See [Preparing data and
configuration](https://mrc-ide.github.io/anatembea/articles/preparing-an-analysis.md)
and [Tuning an anatembea
pMCMC](https://mrc-ide.github.io/anatembea/articles/tuning-pmcmc.md)
first.

The workflow has two gates:

1.  short independent chains confirm the tuned configuration; and
2.  longer independent chains provide the retained scientific inference.

Completing a run is not the same as demonstrating convergence.

## 1. Confirm with independent chains

Start with two or three chains of about 1,000 iterations at the intended
production particle count. Give every chain a distinct seed and, where
possible, a separate scheduler job. Running separate calls also keeps
outputs, logs, and failure states unambiguous.

``` r

confirmation <- data.frame(
  chain = 1:3,
  seed = c(4101L, 4102L, 4103L),
  particles = 200L,
  steps = 1000L
)

confirmation
#>   chain seed particles steps
#> 1     1 4101       200  1000
#> 2     2 4102       200  1000
#> 3     3 4103       200  1000
```

The following code assumes that `production_covariance` was frozen
during tuning. It is not evaluated during the vignette build.

``` r

confirmation_runs <- lapply(seq_len(nrow(confirmation)), function(i) {
  anatembea::run_pmcmc(
    data_raw = analysis_data,
    comparison = "u5",
    initial = "fitted",
    check_flexibility = TRUE,
    proposal_matrix = production_covariance,
    n_particles = confirmation$particles[i],
    n_steps = confirmation$steps[i],
    seed = confirmation$seed[i],
    output_level = "standard",
    save_state = FALSE,
    save_trajectories = TRUE
  )
})
```

`output_level = "standard"` retains the key fitted trajectories needed
to compare chains without retaining the full diagnostic state. Saving
trajectories costs more memory than the minimal tuning runs, so budget
confirmation jobs accordingly.

## 2. Apply burn-in within each chain

Never concatenate chains and remove burn-in only once. Each chain has
its own initial transient.

``` r

parameter_names <- c("log_init_EIR", "volatility")
burn_fraction <- 0.5

chains <- lapply(confirmation_runs, function(result) {
  first_keep <- floor(nrow(result$pars) * burn_fraction) + 1L
  coda::as.mcmc(
    result$pars[first_keep:nrow(result$pars), parameter_names, drop = FALSE]
  )
})
# mcmc.list() accepts a single list whose elements are mcmc objects.
chains <- coda::mcmc.list(chains)
```

Review each chain before combining summaries:

``` r

plot(chains)                         # traces and marginal densities
coda::effectiveSize(chains)         # combined effective sample size
coda::gelman.diag(chains, autoburnin = FALSE)
coda::gelman.plot(chains, autoburnin = FALSE)
```

[`coda::gelman.diag()`](https://rdrr.io/pkg/coda/man/gelman.diag.html)
reports the classical potential scale reduction factor. Where your
analysis stack supports it, also calculate rank-normalized R-hat and
bulk and tail ESS; the [posterior
package](https://mc-stan.org/posterior/) provides these modern
diagnostics. No single threshold should replace trace inspection and
fitted trajectory checks.

Also calculate chain-specific acceptance and ESS so that one weak chain
is not hidden by a combined summary:

``` r

chain_diagnostics <- do.call(rbind, lapply(seq_along(chains), function(i) {
  acceptance <- 1 - coda::rejectionRate(chains[[i]])
  ess <- coda::effectiveSize(chains[[i]])
  data.frame(
    chain = i,
    parameter = names(ess),
    acceptance = acceptance[names(ess)],
    ess = as.numeric(ess)
  )
}))

chain_diagnostics
```

Do not advance if chains occupy different regions, remain strongly
autocorrelated, concentrate unexpectedly at a parameter boundary, or
produce materially different fitted trajectories. Add length or a third
chain only when that action addresses a specific diagnostic concern.

## 3. Freeze the production design

After the confirmation gate, freeze the proposal covariance, its scale,
the particle count, model settings, and output policy. A reasonable
initial design is two to four independent chains of 5,000–10,000
iterations, but confirmation ESS and runtime should determine the actual
length.

Do not adapt the proposal from retained production samples unless the
adapted run is treated as a new stage and independently confirmed again.

``` r

production <- data.frame(
  chain = 1:4,
  seed = c(5101L, 5102L, 5103L, 5104L),
  particles = 200L,
  steps = 5000L
)
```

## 4. Choose what output to retain

Version 1.1 introduces output controls that let the retained object
match the purpose of the run:

| Need | Suggested settings |
|----|----|
| Parameter diagnostics only | `output_level = "minimal"`, no saved trajectories |
| Key prevalence, EIR, emergence, and under-5 incidence paths | `output_level = "standard"`, save trajectories |
| Full internal diagnostic state | `output_level = "diagnostic"`, save trajectories |

`save_state` controls whether the underlying pMCMC state values are
retained; `save_trajectories` controls whether `history` and `times` are
returned. Full diagnostic output can be large. Use it because a
particular diagnostic needs those states, not merely as a default.

For example, a production chain intended for key fitted-trajectory
summaries might use:

``` r

production_run_1 <- anatembea::run_pmcmc(
  data_raw = analysis_data,
  comparison = "u5",
  initial = "fitted",
  check_flexibility = TRUE,
  proposal_matrix = production_covariance,
  n_particles = production$particles[1],
  n_steps = production$steps[1],
  seed = production$seed[1],
  output_level = "standard",
  save_state = FALSE,
  save_trajectories = TRUE
)

stopifnot(
  !is.null(production_run_1$history),
  !is.null(production_run_1$times)
)
```

Keep output settings consistent across chains that will be summarized
together. If only selected chains generate trajectories, document the
selection rule before inspecting their results.

## 5. Parallel execution without losing the gates

Independent datasets, proposal trials, and chains can run concurrently.
Within-run threading should match the scheduler resources actually
requested. Parallel execution shortens elapsed time; it does not change
the dependency order:

1.  validate inputs and execution;
2.  tune particles and covariance;
3.  compare proposal scales;
4.  confirm with short independent chains; and
5.  run production chains.

Avoid multiplying parallelism accidentally. For example, four scheduler
jobs each requesting four threads are already using up to 16 threads;
adding multiple workers inside every job may oversubscribe the
allocation.

## 6. Respond to failures deliberately

| Symptom | First checks | Appropriate response |
|----|----|----|
| Extreme fitted EIR or ODE stalls | Bounds and rejected proposals | Verify bounded `log_init_EIR`; do not simply raise the ODE step limit |
| Particle-filter failures | Parameter values, seed, particle count | Distinguish stochastic particle depletion from deterministic model failure |
| Very low acceptance | Proposal scale and likelihood noise | Reduce scale, improve covariance, or increase particles when justified |
| Very high acceptance with slow movement | Trace range and ESS/hour | Increase proposal scale |
| Singular pilot covariance | Number and movement of retained draws | Extend or repair the pilot; regularize only small numerical defects |
| One unusually slow dataset | Dates, gaps, prevalence, initialization, logs | Diagnose that dataset separately |
| One failed chain | Its seed, logs, and exit state | Retain the failure record; do not silently replace it |

Changing `ode_max_steps`, `ode_atol`, or `ode_rtol` changes numerical
behavior. Treat such changes as part of the recorded analysis
specification and justify them with diagnostics.

## 7. Keep a reconstructible run record

Save raw outputs and a metadata record for every run, including failed
runs.

``` r

production_arguments <- list(
  comparison = "u5",
  initial = "fitted",
  check_flexibility = TRUE,
  proposal_matrix = production_covariance,
  n_particles = production$particles[1],
  n_steps = production$steps[1],
  n_threads = 4,
  n_chains = 1,
  n_workers = 1,
  seed = production$seed[1],
  prop_treated = 0.4,
  seasonality_on = TRUE,
  output_level = "standard",
  save_state = FALSE,
  save_trajectories = TRUE,
  ode_atol = 1e-6,
  ode_rtol = 1e-6,
  ode_max_steps = 1e7
)

run_with_record <- function(data, arguments, input_path = NULL) {
  started_at <- Sys.time()
  captured_warnings <- character()
  error_message <- NULL

  result <- withCallingHandlers(
    tryCatch(
      do.call(anatembea::run_pmcmc, c(list(data_raw = data), arguments)),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  finished_at <- Sys.time()
  exit_state <- if (is.null(error_message)) "completed" else "failed"

  record <- list(
    analysis_id = "tanzania-u5-example",
    input_hash = if (is.null(input_path)) NA_character_ else
      unname(tools::md5sum(input_path)),
    git_commit = tryCatch(
      system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      error = function(e) NA_character_
    ),
    package_versions = vapply(
      c("anatembea", "mcstate", "dust", "odin.dust", "odin"),
      function(x) as.character(utils::packageVersion(x)),
      character(1)
    ),
    run_arguments = arguments,
    scheduler = list(
      job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_)
    ),
    result = list(
      started_at = started_at,
      finished_at = finished_at,
      runtime = if (is.null(result)) NA else result$run_time,
      exit_state = exit_state,
      error = error_message,
      warnings = captured_warnings
    )
  )

  list(result = result, record = record)
}

recorded_run <- run_with_record(analysis_data, production_arguments)
saveRDS(recorded_run$record, "chain-1-run-record.rds")

if (identical(recorded_run$record$result$exit_state, "failed")) {
  stop(recorded_run$record$result$error)
}
production_run_1 <- recorded_run$result
```

Also retain warnings, errors, raw samples, trajectories, diagnostics,
reviewer decisions, submission time, and the rule used to advance or
reject the run. A production result should be reconstructible without an
interactive session or an unversioned package installation.

## Final gate

Report production results only after:

- every chain has been reviewed independently;
- between-chain convergence and within-chain mixing are acceptable;
- bulk and tail behavior are adequate for the reported quantities;
- fitted trajectories are scientifically plausible;
- boundary behavior and failures have been investigated; and
- inputs, code, settings, seeds, outputs, and decisions are archived.

When these conditions are met, tuning has done its job: the final
inference is based on a configuration selected before the retained
production samples were interpreted.

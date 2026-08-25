# Manual benchmark for pMCMC output modes.
#
# This script is intentionally not part of package checks. Run it from the
# package root after installing development dependencies.

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".")
} else {
  library(anatembea)
}

make_benchmark_data <- function() {
  dplyr::filter(
    dplyr::select(anatembea::sim_data_tanzania, month, positive, tested),
    month <= zoo::as.yearmon("Dec 2019")
  )
}

run_mode <- function(label, save_state, save_trajectories, output_level) {
  gc()
  start <- proc.time()[["elapsed"]]
  result <- anatembea::run_pmcmc(
    data_raw = make_benchmark_data(),
    target_prev = 0.4,
    n_particles = 10,
    n_steps = 5,
    n_threads = 1,
    n_chains = 1,
    n_workers = 1,
    save_state = save_state,
    save_trajectories = save_trajectories,
    output_level = output_level
  )
  runtime <- proc.time()[["elapsed"]] - start

  list(
    label = label,
    runtime = runtime,
    object_size = object.size(result),
    mcmc_columns = colnames(result$mcmc),
    pars_columns = colnames(result$pars),
    has_history = !is.null(result$history),
    has_times = !is.null(result$times),
    result = result
  )
}

modes <- list(
  diagnostic = run_mode(
    "diagnostic",
    save_state = TRUE,
    save_trajectories = TRUE,
    output_level = "diagnostic"
  ),
  minimal = run_mode(
    "minimal",
    save_state = FALSE,
    save_trajectories = FALSE,
    output_level = "minimal"
  ),
  standard = run_mode(
    "standard",
    save_state = TRUE,
    save_trajectories = TRUE,
    output_level = "standard"
  )
)

summary <- data.frame(
  mode = names(modes),
  runtime_seconds = vapply(modes, `[[`, numeric(1), "runtime"),
  object_size_mb = vapply(modes, function(x) {
    as.numeric(x$object_size) / 1024^2
  }, numeric(1)),
  has_history = vapply(modes, `[[`, logical(1), "has_history"),
  has_times = vapply(modes, `[[`, logical(1), "has_times")
)

posterior_columns_match <- list(
  mcmc = identical(modes$diagnostic$mcmc_columns, modes$minimal$mcmc_columns) &&
    identical(modes$diagnostic$mcmc_columns, modes$standard$mcmc_columns),
  pars = identical(modes$diagnostic$pars_columns, modes$minimal$pars_columns) &&
    identical(modes$diagnostic$pars_columns, modes$standard$pars_columns)
)

print(summary)
print(posterior_columns_match)

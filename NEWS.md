# anatembea 1.1.0

- Added opt-in pMCMC output modes that reduce runtime and memory use when full
  state trajectories are not required. Existing calls retain the previous
  full-output behaviour by default; use `output_level = "minimal"`,
  `save_state = FALSE`, and `save_trajectories = FALSE` for the smallest output.
- Bounded fitted `log_init_EIR` proposals and corrected their prior
  specification, improving the reliability of `initial = "fitted"` runs.
- Added staged articles for preparing data and configuration, tuning pMCMC,
  and assessing confirmation and production chains.
- Added user controls for ODE absolute tolerance, relative tolerance, and the
  maximum number of ODE steps.

# anatembea 1.0 (accepted manuscript freeze)

- Frozen reproducibility point for the accepted ANC pMCMC manuscript, published in **The Lancet Microbe**.
- Published article: <https://www.sciencedirect.com/science/article/pii/S2666524726000704>.
- Associated preservation tag: `v1.0-paper`.
- Intended manuscript runtime context: R 4.4 with `odin` 1.5.11, `dust` 0.15.3, and `mcstate` 0.9.22.
- Analysis repository: <https://github.com/jt-hicks/anc_pmcmc_pub>.

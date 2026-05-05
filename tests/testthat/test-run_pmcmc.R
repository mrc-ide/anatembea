test_that("data_process formats data correctly", {
  # Create minimal test dataset
  test_data <- data.frame(
    month = zoo::as.yearmon(c("Jan 2015", "Feb 2015", "Mar 2015")),
    positive = c(10, 15, 20),
    tested = c(100, 100, 100)
  )

  result <- data_process(test_data, start_pf_time = 30, check_flexibility = FALSE)

  expect_named(result, c("data", "stochastic_schedule", "start_stoch", "time_origin"))
  expect_true(is.data.frame(result$data) || inherits(result$data, "particle_filter_data"))
  expect_true(inherits(result$time_origin, "Date"))
  expect_true(is.integer(result$stochastic_schedule) || is.numeric(result$stochastic_schedule))
  # time_origin should be Jan 1 of the year before first observation
  expect_equal(format(result$time_origin, "%Y-%m-%d"), "2014-01-01")
})

test_that("format_na handles missing data correctly", {
  df <- data.frame(
    month = zoo::as.yearmon(c("Jan 2015", "Feb 2015")),
    positive = c(10, NA),
    tested = c(100, NA)
  )

  result <- format_na(df)

  expect_equal(result$tested[2], 0)
  expect_equal(result$positive[2], 0)
  expect_equal(result$positive[1], 10)
})

test_that("compare_u5 returns numeric log-likelihood", {
  # state is a matrix (rows = state vars, cols = particles)
  state <- matrix(c(0.3, 0.25, 0.35), nrow = 1)
  rownames(state) <- "prev"
  observed <- list(positive = 30, tested = 100)

  ll <- compare_u5(state, observed)

  expect_true(is.numeric(ll))
  expect_equal(length(ll), 3)
  expect_true(all(ll <= 0))  # log-likelihoods are <= 0 for binomial

  # Missing data should return zeros
  observed_na <- list(positive = NA, tested = 100)
  ll_na <- compare_u5(state, observed_na)
  expect_equal(ll_na, rep(0, 3))
})

test_that("get_odds_from_prev and get_prev_from_log_odds are inverse functions", {
  prev_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  for (p in prev_vals) {
    odds <- get_odds_from_prev(p)
    log_odds <- log(odds)
    p_recovered <- get_prev_from_log_odds(log_odds)
    expect_equal(p_recovered, p, tolerance = 1e-10)
  }
})

test_that("admin_match returns 0 when no country/admin specified", {
  result <- admin_match(admin_unit = NULL, country = NULL,
                        admin_units_seasonal = admin_units_seasonal)
  expect_equal(result, 0)
})

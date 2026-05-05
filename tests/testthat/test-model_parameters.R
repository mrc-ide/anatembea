test_that("Load correlation coefficients", {
  expect_equal(nrow(pgsgmg_corr_sample), 1000)
  expect_equal(pgsgmg_corr_sample[1,'av_lo_child'],-2.369)
})

test_that("model_param_list_create returns a list with expected parameters", {
  mpl <- model_param_list_create(comparison = 'u5', avg_prev = 0.3)

  # Check it's a list
  expect_true(is.list(mpl))

  # Check key epidemiological parameters are present
  expect_true("eta" %in% names(mpl))
  expect_true("b0" %in% names(mpl))
  expect_true("phi0" %in% names(mpl))
  expect_true("DY" %in% names(mpl))

  # Check derived parameters
  expect_true("fv0" %in% names(mpl))
  expect_true("av0" %in% names(mpl))
  expect_true("Surv0" %in% names(mpl))

  # Check log_OR parameters are computed
  expect_true("log_OR_pg_v_c" %in% names(mpl))
  expect_true("log_OR_pm_v_pp" %in% names(mpl))
})

test_that("model_param_list_create errors on duplicate extra params", {
  # Named params like eta can be overridden without error
  expect_error(
    model_param_list_create(comparison = 'u5', avg_prev = 0.3, eta = 0.0002),
    NA  # NA means expect no error
  )

  # Extra param with same name as a derived list element should error
  expect_error(
    model_param_list_create(comparison = 'u5', avg_prev = 0.3, fv0 = 999)
  )
})

test_that("model_param_list_create handles flex comparison", {
  mpl <- model_param_list_create(comparison = '2to10', avg_prev = 0.3,
                                 target_prev = 0.3)

  expect_equal(mpl$age_min, 2)
  expect_equal(mpl$age_max, 10)
})

test_that("get_anc_from_u5 returns prevalence values in [0,1]", {
  result <- get_anc_from_u5(prev_u5 = 0.3, avg_prev = 0.2)

  expect_true(is.list(result))
  expect_named(result, c("prev_preg_pg", "prev_preg_sg", "prev_preg_mg", "prev_preg_all"))
  for (nm in names(result)) {
    expect_true(result[[nm]] >= 0 && result[[nm]] <= 1,
                info = paste(nm, "should be in [0,1]"))
  }
  # Multigravidae should have lower prevalence than primigravidae
  expect_lt(result$prev_preg_mg, result$prev_preg_pg)
})


test_that("Load correlation coefficients", {
  corr_sample <- mamasante::load_file("pgsgmg_corr_sample.RDS")
  expect_equal(nrow(corr_sample), 1000)
  expect_equal(corr_sample[1,'av_lo_child'],-2.369)
})

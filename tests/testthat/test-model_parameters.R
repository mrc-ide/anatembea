test_that("Load correlation coefficients", {
  expect_equal(nrow(pgsgmg_corr_sample), 1000)
  expect_equal(pgsgmg_corr_sample[1,'av_lo_child'],-2.369)
})

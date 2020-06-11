test_that("works", {
  fit <- blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4)
  expect_s3_class(fit, "blblm_parallel")
  beta <- coef(fit)
  expect_equal(length(beta), 2)
})

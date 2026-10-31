test_that("hyper.like returns weighted likelihoods and sigma offsets", {
  X <- matrix(c(0, 3, 1, 3), ncol = 2, byrow = TRUE)
  covarray <- array(0, dim = c(2, 2, nrow(X)))
  parm <- c(0, 2, 2)
  point_log_likelihood <- -0.5 * (log(4) + 1 / 4)

  expect_equal(
    hyper.like(parm, X, covarray, weights = c(1, 2), output = "sum"),
    3 * point_log_likelihood
  )
  expect_equal(
    hyper.like(parm, X, covarray, output = "val"),
    rep(point_log_likelihood, 2)
  )
  expect_equal(hyper.like(parm, X, covarray, output = "sig"), c(0.5, 0.5))
})

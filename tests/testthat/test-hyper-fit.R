test_that("hyper.fit returns fitted parameters and predictions", {
  x <- seq(-2, 2, length.out = 10)
  y <- 1.5 * x + 0.7 +
    c(-0.2, 0.1, 0.15, -0.1, 0.05, -0.05, 0.08, -0.12, 0.18, -0.03)
  X <- cbind(x = x, y = y)

  fit <- hyper.fit(X, parm = c(1.5, 0.7, 0.5))

  expect_s3_class(fit, "hyper.fit")
  expect_equal(fit$N, nrow(X))
  expect_equal(fit$dims, ncol(X))
  expect_length(fit$predict.vert.axis, nrow(X))
  expect_true(is.finite(fit$LL$sum))
})

test_that("hyper.fit rejects data with fewer than two dimensions", {
  expect_error(
    hyper.fit(matrix(seq_len(4), ncol = 1)),
    "2 or more columns"
  )
})

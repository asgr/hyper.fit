test_that("rotation matrices rotate vectors by the requested angle", {
  expect_equal(
    drop(makerotmat2d(90) %*% c(1, 0)),
    c(0, 1),
    tolerance = 1e-12
  )
  expect_equal(
    drop(makerotmat3d(90, dim = "z") %*% c(1, 0, 0)),
    c(0, 1, 0),
    tolerance = 1e-12
  )
})

test_that("covariance helpers preserve covariance and projection values", {
  covmat <- makecovmat2d(sx = 2, sy = 3, corxy = 0.5)
  expect_equal(covmat, matrix(c(4, 3, 3, 9), nrow = 2))
  expect_equal(projcovmat(covmat, c(1, 2)), 52)

  covarray <- makecovarray2d(
    sx = c(2, 1),
    sy = c(3, 2),
    corxy = c(0.5, 0)
  )
  expect_equal(dim(covarray), c(2, 2, 2))
  expect_equal(projcovarray(covarray, c(1, 2)), c(52, 17))
})

test_that("array-vector multiplication checks conformable dimensions", {
  array <- array(seq_len(12), dim = c(2, 3, 2))
  expect_equal(
    arrayvecmult(array, c(1, 2, 3)),
    rbind(c(22, 28), c(58, 64))
  )
  expect_error(
    arrayvecmult(array, c(1, 2)),
    "Second dimension of array must be the same length as vec"
  )
})

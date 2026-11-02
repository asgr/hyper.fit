test_that("orthogonal scatter conversion retains its name and value", {
  original <- hyper.convert(coord = 2, beta = 1, scat = 4)
  normal <- convert(original, coord.type = "normvec", scat.type = "orth")

  expect_equal(normal$scat, c(scat.orth = 4 / sqrt(5)))
})

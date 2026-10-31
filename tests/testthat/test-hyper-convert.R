test_that("hyperplane parameters round-trip between coordinate systems", {
  original <- hyper.convert(coord = 2, beta = 1, scat = 4)
  expect_s3_class(original, "hyper.plane.param")

  normal <- convert(original, coord.type = "normvec", scat.type = "orth")
  expect_equal(normal$coord, c(normvec1 = -0.4, normvec2 = 0.2))
  expect_equal(normal$scat, 4 / sqrt(5))

  round_trip <- convert(
    normal,
    coord.type = "alpha",
    scat.type = "vert.axis",
    vert.axis = 2
  )
  expect_equal(round_trip$coord, original$coord)
  expect_equal(round_trip$beta, original$beta)
  expect_equal(round_trip$scat, original$scat)
})

test_that("theta coordinates represent the corresponding gradients", {
  converted <- hyper.convert(
    coord = 2,
    beta = 1,
    scat = 0.5,
    out.coord.type = "theta"
  )
  expect_equal(unname(converted$coord), atan(2) * 180 / pi)
})

test_that("hyper.convert requires a parameter or coordinate input", {
  expect_error(hyper.convert(), "Must provide either parm or coord input")
})

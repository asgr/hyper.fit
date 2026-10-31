# Agent guidance

## Package overview

`hyper.fit` is an R package for fitting N-dimensional hyperplanes while accounting
for measurement errors and intrinsic scatter.

- Package functions are in `R/`; their help pages are in `man/`.
- Example datasets are in `data/`, and package metadata and dependencies are in
  `DESCRIPTION` and `NAMESPACE`.
- `hyper.convert()` changes hyperplane parameterisations, `hyper.like()`
  evaluates the likelihood, and `hyper.fit()` performs a fit.
- The package also provides geometry, covariance, and plotting helpers.

## Development and tests

- Add unit tests under `tests/testthat/`; `tests/testthat.R` is the package test
  entry point.
- `testthat` is a suggested dependency. Declare additional package dependencies
  in `DESCRIPTION`, not in ad hoc setup scripts.
- Run the tests with `R -q -e 'testthat::test_local()'`.
- Run a package check with `R CMD build .` followed by
  `R CMD check --no-manual hyper.fit_*.tar.gz`.
- Keep tests deterministic. Set a seed when testing randomized helpers, and
  avoid opening graphics devices for tests unless plotting behavior is the
  subject of the test.

## CI

`.github/workflows/main.yml` runs `R CMD check` across operating systems and R
versions. CI uses `r-lib/actions/setup-r-dependencies`; make sure test and
documentation dependencies are available through package metadata or the
workflow's `extra-packages` setting.

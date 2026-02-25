# hyper.fit (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/hyper.fit/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Core package containing all the tools for N-dimensional hyperplane fitting with errors. In simple terms this allows the user to produce robust 1D linear fits for 2D x vs y type data, and robust 2D plane fits to 3D x vs y vs z type data. This hyperplane fitting works generically for any N-1 hyperplane model being fit to a N dimension dataset. All fits include intrinsic scatter in the generative model orthogonal to the hyperplane.

## Installation

### Getting hyper.fit

Currently **hyper.fit** is on CRAN, which means you can install easily with:

```R
install.packages('hyper.fit')
```

If that does not work for some reason, source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/hyper.fit")
library(hyper.fit)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools (which I do not use as the default since it has a few more dependencies, and is tricky to install on HPCs):

```R
install.packages('devtools')
devtools::install_github("asgr/hyper.fit")
library(hyper.fit)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("asgr/hyper.fit")
```

I also have these options set by default in my .Rprofile, which seems to help with some of the remote install issues some people face:

```R
options(download.file.method = "libcurl")
options(repos="http://cran.rstudio.com/")
options(rpubs.upload.method = "internal")
```

If all of these do not work than the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install hyper.fit_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

If none of the above works then you should consider burning your computer in sacrifice to the IO Gods. Then buy a newer *better* computer, and try all the above steps again.

Failing all of the above, please email me for help (or perhaps raise an Issue here, if it really does not seem like a local issue).

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **hyper.fit**:

```R
install.packages(c('magicaxis', 'MASS', 'rgl', 'LaplacesDemon')) # Required packages
install.packages('remotes')
remotes::install_github("asgr/hyper.fit")
```

Assuming you have installed all of the packages that you need/want, you should now be able to load **hyper.fit** within **R** with the usual:

```R
library(hyper.fit)
```

## Code Example

```R
# A very simple 2D example: fit a line to x vs y data

simpledata = cbind(x=1:10, y=c(12.2, 14.2, 15.9, 18.0, 20.1, 22.1, 23.9, 26.0, 27.9, 30.1))
simpfit = hyper.fit(simpledata)
summary(simpfit)
plot(simpfit)
```

To find more long-form examples, including complicated fitting use-cases, please check the help documentation. You can browse these with:

```R
?hyper.fit
```

## Motivation

A general purpose hyperplane fitting package intended to be used for 1D line fitting (2D data) or 2D plane fitting (3D data) and beyond. It was originally developed to support fitting of galaxy scaling relations in astronomy, but it serves as a general purpose regression tool with intrinsic scatter in its own right.

## Contributors

Code:

Aaron Robotham

## References

Main introduction:

Robotham A.S.G., & Obreschkow D., 2015, PASA, 32, 33

## Resources

<https://ui.adsabs.harvard.edu/abs/2015PASA...32...33R/abstract>

## Forums

Please sign up to <http://profit.freeforums.net/> if you want to ask a question (or browse the questions asked).

## License

GPL-3

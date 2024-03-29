##
## RUnit.r - run unit tests in inst/unittests
##
## This file is based on code in the CRAN package 'fBasics', which in
## turn is based on code from the 'gdata' package.
##

pkg <- "cmaes"

if (require("RUnit", quietly = TRUE)) {
  wd <- getwd()
  library(package = pkg, character.only = TRUE)
  path <- system.file("unittests", package = pkg)
  stopifnot(file.exists(path), file.info(path.expand(path))$isdir)
  source(file.path(path, "runner.r"), echo = TRUE)
}

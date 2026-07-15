# merDeriv 0.2-5

## Bug fixes

* `vcov.lmerMod(full = TRUE)` and `estfun.lmerMod()` returned the random
  effect (co)variance results in an order that did not match their labels
  when a random effect block was 3x3 or larger, so results depended on the
  order in which random effects were entered in the model formula
  (Github issue #2, reported by James Pustejovsky). The parameter index
  matrix is now built from `Lambdat`, which `Lind` actually indexes,
  instead of `Lambda`.

* `estfun.glmerMod()` and `llcont.glmerMod()` failed with
  "subscript out of bounds" when the cluster id variable was not a
  consecutive integer sequence starting at 1 (e.g., character ids or
  arbitrary numeric ids; Github issue #3, reported by teindor). Cluster
  rows are now looked up by name rather than by coerced numeric index.

* `vcov.lmerMod(full = TRUE, information = "observed")` failed for REML
  fits because the internal Hessian entries were not coerced to numeric.

* Internal `class(x) == "..."` comparisons were replaced with
  `inherits()`, avoiding "condition has length > 1" errors on
  R >= 4.2.0 for models whose family was specified as a function.

## Other changes

* Test suite migrated to tinytest, with new regression tests covering
  Github issues #2 and #3, and updated lavaan model syntax
  (Github issue #4, reported by Yves Rosseel).

* Documentation fixes, including the `bread.glmerMod()` man page
  (contributed by Ben Bolker).

* Suggested packages (mirt, lmeInfo, nlme, smdata) are now used
  conditionally in tests, and the duplicated Suggests field in
  DESCRIPTION was merged.

# merDeriv 0.2-4

* Released to CRAN 2022-06-23. `estfun.glmerMod()` support for
  binomial and poisson families fit with `glmer()`, with scores of
  cluster-wise log-likelihoods computed via adaptive Gauss-Hermite
  quadrature.

library("lme4")

## Github issue #2: with 3 or more random effects per grouping factor,
## the variance component block of vcov(., full = TRUE) and the
## variance component columns of estfun() came back in an order that
## did not match their labels, so results depended on the order in
## which random effects were entered in the formula.
set.seed(42)
J <- 60; nj <- 10; N <- J * nj
G <- matrix(c(4, 1, .8,  1, 2, .5,  .8, .5, 1.5), 3, 3)
Gchol <- chol(G)
simd <- data.frame(g = rep(1:J, each = nj), x1 = rnorm(N), x2 = rnorm(N))
re <- matrix(rnorm(J * 3), J, 3) %*% Gchol
simd$y <- 2 + .5 * simd$x1 - .3 * simd$x2 + re[simd$g, 1] +
  re[simd$g, 2] * simd$x1 + re[simd$g, 3] * simd$x2 + rnorm(N, 0, 1.5)

fit1 <- lmer(y ~ x1 + x2 + (x1 + x2 | g), simd, REML = FALSE)
fit2 <- lmer(y ~ x1 + x2 + (x2 + x1 | g), simd, REML = FALSE)

## map parameter names of fit2 onto fit1 (covariance labels may have
## their two variable names in either order, so compare sorted labels)
normnames <- function(v) {
  ifelse(grepl("^cov_g\\.", v),
         paste0("VC:", sapply(strsplit(sub("^cov_g\\.", "", v), "\\.",
                                       fixed = FALSE),
                              function(s) paste(sort(s), collapse = "+"))),
         v)
}

for (info in c("expected", "observed")) {
  V1 <- as.matrix(vcov(fit1, full = TRUE, information = info))
  V2 <- as.matrix(vcov(fit2, full = TRUE, information = info))
  idx <- match(normnames(colnames(V1)), normnames(colnames(V2)))
  expect_false(any(is.na(idx)))
  expect_equal(V1, V2[idx, idx], check.attributes = FALSE, tolerance = 1e-3,
               info = paste("vcov consistency across RE orderings,", info))
}

S1 <- estfun(fit1)
S2 <- estfun(fit2)
idx <- match(normnames(colnames(S1)), normnames(colnames(S2)))
expect_equal(unname(as.matrix(S1)), unname(as.matrix(S2))[, idx],
             check.attributes = FALSE, tolerance = 1e-3,
             info = "estfun consistency across RE orderings")

## the variance component block must also be correctly labeled: the
## diagonal variance entries follow theta order (column-major lower
## triangle), so the entry labeled cov_g.x1 must be the sampling
## variance of the x1 random effect variance, which is invariant to
## formula order once matched by label.
v1 <- diag(as.matrix(vcov(fit1, full = TRUE)))
v2 <- diag(as.matrix(vcov(fit2, full = TRUE)))
expect_equal(v1["cov_g.x1"], v2["cov_g.x1"], tolerance = 1e-3,
             check.attributes = FALSE)
expect_equal(v1["cov_g.x2"], v2["cov_g.x2"], tolerance = 1e-3,
             check.attributes = FALSE)

## Github issue #3: estfun.glmerMod() failed with "subscript out of
## bounds" when cluster ids were numeric but not consecutive integers
## starting at 1 (and for cluster ids that are not numbers at all).
set.seed(7)
ng <- 40
gd <- data.frame(cID = rep(sample(1000:9999, ng), each = 3),
                 x = rnorm(3 * ng))
gd$y <- rbinom(nrow(gd), 1, plogis(.5 + .4 * gd$x + rep(rnorm(ng), each = 3)))

gfit <- glmer(y ~ x + (1 | cID), data = gd, family = binomial, nAGQ = 5L)
gs <- estfun(gfit)
expect_true(is.matrix(gs))
expect_equal(nrow(gs), ng)
expect_true(all(abs(colSums(gs)) < .05))

## character cluster ids
gd$cIDc <- paste0("id", gd$cID)
gfitc <- glmer(y ~ x + (1 | cIDc), data = gd, family = binomial, nAGQ = 5L)
gsc <- estfun(gfitc)
expect_equal(unname(gs), unname(gsc), tolerance = 1e-6)

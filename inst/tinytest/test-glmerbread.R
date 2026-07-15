library(merDeriv)
if (require("smdata")) {
  data(finance, package = "smdata")

   lme4fit <- lme4::glmer(corr ~ jmeth + (1 | item), data = finance,
                 family = binomial, nAGQ = 20)

  # bread component for all parameters
  b <- sandwich::bread(lme4fit, full = TRUE, ranpar = "var")
  print(class(b))
  stopifnot(inherits(b,"matrix"),
            identical(colnames(b),
         c("(Intercept)", "jmeth2ci", "jmeth3ei", "cov_item.(Intercept)")))
}


library("lme4")
library("mirt")
library("reshape2")

data(VerbAgg)

## check llcont.glmerMod:
fmVA0 <- glmer(r2 ~ Anger + (Gender | item), family = binomial, data = VerbAgg, nAGQ=0L)
ll <- llcont.glmerMod(fmVA0)
expect_equal(as.numeric(round(sum(ll),4)), as.numeric(round(logLik(fmVA0),4)))

## with different nAGQ
ll <- llcont.glmerMod(fmVA0, nAGQ=5)
expect_equal(as.numeric(round(sum(ll))), as.numeric(round(logLik(fmVA0)))) # (these disagree a little)

## error for multiple clustering variables
fm2 <- glmer(r2 ~ Anger + (Anger | item) + (1 | id), family = binomial, data = VerbAgg, nAGQ=0L)
expect_error(llcont.glmerMod(fm2))

## compare with score produced by mirt
## score from lme4
data <- expand.table(LSAT7)
data$person <- as.numeric(rownames(data))
datalong <- melt(data, id = c("person"))
fit <- glmer(value ~ -1 + variable + (1|person), family = binomial, 
             data = datalong, nAGQ = 30)
score <- estfun.glmerMod(fit)

# score from mirt
data <- expand.table(LSAT7)
mod <- mirt(data, 1, itemtype="Rasch")
sc1 <- estfun.AllModelClass(mod)
## compare
expect_true(all(abs(score[,1:5] - sc1[,1:5]) < .01))


## non-canonical link check
fit2 <- glmer(value ~ -1 + variable + (1|person),
              family = binomial(link = "cloglog"), 
              data = datalong, nAGQ = 5)
fixscore2 <- estfun.glmerMod(fit2)
expect_true(all(abs(colSums(fixscore2)) < .1))

## check poisson 
## fit poisson
set.seed(1090)
simfun <- function(ng = 20, nr = 100, fsd = 1, indsd = 0.2, b = c(1,2)) {
  ntot <- nr * ng
  b.reff <- rnorm(ng, sd = fsd)
  b.rind <- rnorm(ntot, sd = indsd)
  x <- runif(ntot)
  dd <- data.frame(x, f = factor(rep(c(1:ng), each = nr)),
                   obs = 1:ntot)
  dd$eta0 <- model.matrix(~x, data = dd) %*% b
  dd$eta <- with(dd, eta0 + b.reff[f] + b.rind[obs])
  dd$mu <- exp(dd$eta)
  dd$y <- with(dd, rpois(ntot, lambda = mu))
  dd
}
dd <- simfun()
m0 <- glmer(y ~ x + (1|f), family = "poisson", data = dd, nAGQ = 20)
score <- estfun.glmerMod(m0)
expect_true(all(abs(colSums(score)) < .3))




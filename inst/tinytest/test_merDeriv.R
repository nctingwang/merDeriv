
library("lme4")
library("mirt")
library("lavaan")

## lmm models checks.
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

# log likelihood check
expect_equal(as.numeric(round(sum(llcont(lme4fit, level = 1)),4)), 
             as.numeric(round(logLik(lme4fit),4)))

## score check
score2 <- estfun(lme4fit, ranpar = "var", level = 2)

testwide <- as.data.frame(t(matrix(sleepstudy$Reaction, 10, length(unique(sleepstudy$Subject)))))
names(testwide) <- paste("d", 1:10, sep = "")
testwide$Subject <- unique(sleepstudy$Subject)

## describe latent model
latent <- 'i =~ 1*d1 + 1*d2 + 1*d3 + 1*d4 + 1*d5
+ 1*d6 + 1*d7 + 1*d8 + 1*d9 + 1*d10

s = ~ 0*d1 + 1*d2 + 2*d3 + 3*d4 + 4*d5
+ 5*d6 + 6*d7 + 7*d8 + 8*d9 + 9*d10

d1 ~~ evar*d1
d2 ~~ evar*d2
d3 ~~ evar*d3
d4 ~~ evar*d4
d5 ~~ evar*d5
d6 ~~ evar*d6
d7 ~~ evar*d7
d8 ~~ evar*d8
d9 ~~ evar*d9
d10 ~~ evar*d10

## reparameterize as sd
sdevar := sqrt(evar)
i ~~ ivar*i
isd := sqrt(ivar)
s ~~ svar*s
ssd := sqrt(svar)
i ~~ iscov*s
rho := iscov/(isd*ssd)'
## fit model in lavaan
lavaanfit <- growth(latent, data = testwide, estimator = "ML")
score2_lavaan <- as.data.frame(lavaan::estfun.lavaan(lavaanfit))
score2_lavaan <- score2_lavaan[,c(5, 6, 2, 4, 3, 1)]
expect_true(all(score2_lavaan-score2 < 0.001))

## vcov and bread dimension check
expect_equal(dim(vcov(lme4fit)), c(2,2))
expect_equal(dim(vcov(lme4fit, full = TRUE, information = "observed")), c(6,6))
expect_equal(dim(bread(lme4fit, full = TRUE)), c(6,6))

## glmm models checks
## check llcont.glmerMod:
fmVA0 <- glmer(r2 ~ Anger + (Gender | item), family = binomial, data = VerbAgg, nAGQ=0L)
ll <- llcont.glmerMod(fmVA0)
expect_equal(as.numeric(round(sum(ll),4)), as.numeric(round(logLik(fmVA0),4)))

## with different nAGQ
ll <- llcont.glmerMod(fmVA0, nAGQ=5)
expect_equal(as.numeric(round(sum(ll))), as.numeric(round(logLik(fmVA0)))) # (these disagree a little)

## error for multiple clustering variable
data(VerbAgg)
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




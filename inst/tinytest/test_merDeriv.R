library("lme4")

## lmm models checks.
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

# log likelihood check
expect_equal(as.numeric(round(sum(llcont(lme4fit, level = 1)),4)), 
             as.numeric(round(logLik(lme4fit),4)))

## score check
score2 <- estfun(lme4fit, ranpar = "var", level = 2)
expect_true(is.matrix(score2) && nrow(score2) == 18L)

## vcov check
expect_true(dim(vcov(lme4fit, full = TRUE))[1] == 6L)

## cluster of size 1
sub1 <- which(sleepstudy$Subject == unique(sleepstudy$Subject)[1])
ss2 <- sleepstudy[-(sub1[-1]),]
subfit <- lmer(Reaction ~ Days + (Days|Subject), data=ss2, REML=FALSE)

expect_equal(as.numeric(round(sum(llcont(subfit, level = 1)),4)), 
             as.numeric(round(logLik(subfit),4)))

expect_equal(class(estfun(subfit))[1], "matrix")

## vcov and bread dimension check
expect_equal(dim(vcov(lme4fit)), c(2,2))
expect_equal(dim(vcov(lme4fit, full = TRUE, information = "observed")), c(6,6))
expect_equal(dim(bread(lme4fit, full = TRUE)), c(6,6))

## glmm models checks
fmVA0 <- glmer(r2 ~ Anger + (Gender | item), family = binomial, data = VerbAgg, nAGQ=0L)
expect_error(vcov(fmVA0, full = TRUE))

fmVA0 <- glmer(r2 ~ Anger + (Gender | item), family = binomial, data = VerbAgg, nAGQ=1L)
ll <- llcont.glmerMod(fmVA0)
expect_equal(as.numeric(round(sum(ll),4)), as.numeric(round(logLik(fmVA0),4)))
expect_true(dim(vcov(fmVA0, full = TRUE))[1] == 5L)

## with different nAGQ
ll <- llcont.glmerMod(fmVA0, nAGQ=5)
expect_equal(as.numeric(round(sum(ll))), as.numeric(round(logLik(fmVA0)))) # (these disagree a little)

## error for multiple clustering variable
data(VerbAgg)
fm2 <- glmer(r2 ~ Anger + (Anger | item) + (1 | id), family = binomial, data = VerbAgg, nAGQ=0L)
expect_error(llcont.glmerMod(fm2))
expect_error(vcov(fm2, full=TRUE))

## compare with score produced by mirt
## score from lme4
if (requireNamespace("mirt", quietly = TRUE)) {
library("mirt")
data <- expand.table(LSAT7)
data$person <- as.numeric(rownames(data))
datalong <- as.data.frame(matrix(NA, nrow(data)*5, 3))
colnames(datalong) <- c("person", "variable", "value")
datalong$variable <- rep(paste0("Item.", seq(1,5)), 
                         each = nrow(data))
datalong$person <- rep(as.numeric(rownames(data)),5)
datalong$value <- as.numeric(as.vector(as.matrix(data[,c(1:5)])))
#datalong <- melt(data, id = c("person"))
fit <- glmer(value ~ -1 + variable + (1|person), family = binomial, 
             data = datalong, nAGQ = 30)
score <- estfun.glmerMod(fit)
expect_true(dim(vcov(fit, full=TRUE))[1] == 6)

# score from mirt
#data <- expand.table(LSAT7)
#mod <- mirt(data, 1, itemtype="Rasch")
#sc1 <- estfun.AllModelClass(mod)
## compare
#expect_true(all(abs(score[,1:5] - sc1[,1:5]) < .01))


## non-canonical link check
fit2 <- glmer(value ~ -1 + variable + (1|person),
              family = binomial(link = "cloglog"), 
              data = datalong, nAGQ = 5)
fixscore2 <- estfun.glmerMod(fit2)
expect_true(all(abs(colSums(fixscore2)) < .1))
expect_true(dim(vcov(fit2, full=TRUE))[1] == 6)
}

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
expect_true(dim(vcov(m0, full=TRUE))[1] == 3)

## check multiple groups vcov
if (requireNamespace("nlme", quietly = TRUE) &&
    requireNamespace("lmeInfo", quietly = TRUE)) {
library("nlme")
library("lmeInfo")
data(Oats, package = "nlme")
Oats.lme <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, method="REML", data=Oats)
Oats.lmer <- lmer(yield ~ nitro*Variety+(1|Block/Variety), REML=TRUE, data=Oats)

nlmeres <- c(varcomp_vcov(Oats.lme)[2,2], varcomp_vcov(Oats.lme)[3,3])
lme4res <- c(vcov.lmerMod(Oats.lmer, full = TRUE, ranpar = "var")[7, 7],
             vcov.lmerMod(Oats.lmer, full = TRUE, ranpar = "var")[9, 9])

expect_true(all(abs(nlmeres-lme4res) < 1))
}

## check multiple groups score
nestmod <- lmer(strength ~ 1 + (1|sample) + (1|batch), Pastes, REML = TRUE)
crossmod <- lmer(diameter ~ 1 + (1|plate) + (1|sample),
                 Penicillin, REML = TRUE)

expect_true(sum(dim(estfun.lmerMod(nestmod, level = "cluster"))==c(40,4))==2)
expect_true(sum(dim(estfun.lmerMod(crossmod, level = "cluster"))==c(30,4))==2)

if (requireNamespace("nlme", quietly = TRUE)) {
data(Oats, package = "nlme")
Oats.lmer <- lmer(yield ~ nitro*Variety+(1|Block/Variety), REML=TRUE, data=Oats)
expect_true(sum(dim(estfun.lmerMod(Oats.lmer, level = "cluster"))==c(24,9))==2)
}


## Gamma example
# set.seed(983)
# dd <- cbind.data.frame(y=rgamma(100,1,2), x1=rnorm(100,5), x2=rnorm(100,3),
#                        id=rep(1:100,each=10))
# m0 <- glmer(y ~ x1 + x2 + (1|id), family=Gamma, data=dd)
# expect_true(abs(sum(llcont.glmerMod(m0)) - logLik(m0)) < .001)

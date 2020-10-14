
library("lme4")
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

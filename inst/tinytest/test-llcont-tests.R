library("lme4")
m1 <- lmer(Reaction ~ 1 + (1|Subject), sleepstudy)
m1_ML <- refitML(m1)
m1_ML_wts <- update(m1_ML, weights = rep(1, nrow(sleepstudy)))
m1_glmer <- glmer(Reaction ~ 1 + (1|Subject), sleepstudy,
                 family = gaussian(link = "log"))

expect_error(llcont(m1), "only works for ML estimation")
expect_error(llcont(m1_ML_wts), "weights specification")
expect_error(llcont(m1_glmer), "family has to be binomial or poisson")

# OR vs RR in R: base R first, packages optional
#
# KEY MESSAGE: epiR (or any package) is NOT required. A fitted glm IS the
# model; exponentiating its coefficients gives OR or RR. Packages only
# add convenience (pretty tables, robust SEs, one-shot 2x2 summaries).

set.seed(42)
n <- 2000
smoker <- rbinom(n, 1, 0.3)
hypertension <- rbinom(n, 1, plogis(-2.5 + 0.8 * smoker))
dat <- data.frame(smoker, hypertension)

# --- 1) base glm: the model IS the difference ---

# logistic regression (logit link) -> odds ratio
fit_or <- glm(hypertension ~ smoker, binomial, dat)
round(exp(cbind(OR = coef(fit_or), confint.default(fit_or))), 3)

# log-binomial (log link) -> risk ratio.
# poisson(link = "log") is the classic trick: Poisson SEs are wrong for
# binary outcomes, so pair it with robust SEs (see section 3).
fit_rr <- glm(hypertension ~ smoker, poisson(link = "log"), dat)
round(exp(cbind(RR = coef(fit_rr), confint.default(fit_rr))), 3)

# --- 2) raw 2x2 table in base R (what epi.2by2 would print) ---
tab <- table(smoker, hypertension)
a <- tab[2, 2]; b <- tab[2, 1]; c <- tab[1, 2]; d <- tab[1, 1]
OR <- (a * d) / (b * c)
RR <- (a / (a + b)) / (c / (c + d))
round(c(OR = OR, RR = RR), 3)   # matches the glm estimates above

# --- 3) OPTIONAL: robust SEs for log-binomial RR ---
# install.packages("sandwich")
library(sandwich)
se <- sqrt(diag(vcovHC(fit_rr)))
round(cbind(RR = exp(coef(fit_rr)),
            lo = exp(coef(fit_rr) - 1.96 * se),
            hi = exp(coef(fit_rr) + 1.96 * se)), 3)



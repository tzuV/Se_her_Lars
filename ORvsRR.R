set.seed(42)
n <- 2000

# --- generate data ---
age    <- rnorm(n, 50, 12)
smoker <- rbinom(n, 1, 0.3)
hypertension <- rbinom(n, 1, plogis(-3 + 0.03 * age + 0.8 * smoker))

etas <- cbind(0, -1.0 + 0.02 * age + 0.5 * smoker,
              -1.5 + 0.04 * age + 0.9 * smoker)
probs <- exp(etas) / rowSums(exp(etas))
severity <- factor(apply(probs, 1, function(p) sample(0:2, 1, prob = p)),
                   labels = c("mild", "moderate", "severe"))

dat <- data.frame(age, smoker, hypertension, severity)

# --- OR: binary ---
fit_bin <- glm(hypertension ~ smoker + age, binomial, dat)
round(exp(cbind(OR = coef(fit_bin), confint.default(fit_bin))), 3)

# --- OR: multinomial ---
library(nnet)
fit_multi <- multinom(severity ~ smoker + age, dat, trace = FALSE)
round(exp(coef(fit_multi)), 3)

# --- RR: binary (log-binomial) ---
fit_rr <- glm(hypertension ~ smoker + age, poisson(link = "log"), dat)
round(exp(cbind(RR = coef(fit_rr), confint.default(fit_rr))), 3)

# --- RR: multinomial (standardized over age) ---
p1 <- colMeans(predict(fit_multi, transform(dat, smoker = 1), "probs"))
p0 <- colMeans(predict(fit_multi, transform(dat, smoker = 0), "probs"))
round(p1 / p0, 3)
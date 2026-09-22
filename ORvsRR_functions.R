# OR and RR built from scratch: the actual formulas
#
# 2x2 table layout (rows = exposure, cols = outcome):
#            disease   no disease
#  exposed     a          b
#  unexposed   c          d
#
#  OR = (a/b) / (c/d) = (a*d)/(b*c)
#     odds of disease in exposed vs odds in unexposed
#
#  RR = [a/(a+b)] / [c/(c+d)]
#     risk (probability) of disease in exposed vs unexposed
#
# Variance (Wald, on the log scale) and 95% CI:
#  Var(log OR) = 1/a + 1/b + 1/c + 1/d
#  Var(log RR) = 1/a - 1/(a+b) + 1/c - 1/(c+d)

or2x2 <- function(tab, conf.level = 0.95) {
  a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
  est <- (a * d) / (b * c)
  se <- sqrt(1/a + 1/b + 1/c + 1/d)
  z <- qnorm(1 - (1 - conf.level) / 2)
  c(OR = est, lower = est * exp(-z * se), upper = est * exp(z * se))
}

rr2x2 <- function(tab, conf.level = 0.95) {
  a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
  r1 <- a / (a + b)   # risk in exposed
  r0 <- c / (c + d)   # risk in unexposed
  est <- r1 / r0
  se <- sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))
  z <- qnorm(1 - (1 - conf.level) / 2)
  c(RR = est, lower = est * exp(-z * se), upper = est * exp(z * se))
}

# --- simple examples ---

# Example 1: classic 2x2, moderate association
tab1 <- matrix(c(30, 70, 15, 85), nrow = 2, byrow = TRUE,
               dimnames = list(exposed = c("yes", "no"), disease = c("yes", "no")))
round(or2x2(tab1), 3)   # OR = (30*85)/(70*15) = 2.429
round(rr2x2(tab1), 3)   # RR = (30/100)/(15/100) = 2.000
# OR > RR here because the outcome is not rare; both point the same direction.

# Example 2: rare outcome -> OR and RR nearly identical
tab2 <- matrix(c(5, 995, 2, 998), nrow = 2, byrow = TRUE,
               dimnames = list(exposed = c("yes", "no"), disease = c("yes", "no")))
round(or2x2(tab2), 3)   # OR  = 2.506
round(rr2x2(tab2), 3)   # RR  = 2.499
# The "rare disease assumption": OR approximates RR when risk < ~10%.

# Example 3: no association -> both equal 1
tab3 <- matrix(c(10, 90, 20, 180), nrow = 2, byrow = TRUE,
               dimnames = list(exposed = c("yes", "no"), disease = c("yes", "no")))
round(or2x2(tab3), 3)   # 1
round(rr2x2(tab3), 3)   # 1

# Example 4: one more dimension -> 2x2x2 (exposure x disease x sex)
#
# Stratify on a third variable and estimate OR/RR per stratum. To get ONE
# number for the whole table, pool the strata with Mantel-Haenszel:
#
#   MH OR = sum(a_i * d_i / n_i) / sum(b_i * c_i / n_i)
#   MH RR = sum(a_i * (c_i + d_i) / n_i) / sum(c_i * (a_i + b_i) / n_i)
#
# where i indexes strata and n_i is the stratum size.

or_mh <- function(tab3d) {
  w <- function(i) {                              # cell counts of stratum i
    a <- tab3d[i, 1, 1]; b <- tab3d[i, 1, 2]
    c <- tab3d[i, 2, 1]; d <- tab3d[i, 2, 2]
    c(a, b, c, d, n = a + b + c + d)
  }
  k <- dim(tab3d)[1]
  cells <- t(vapply(1:k, w, numeric(5)))
  num <- sum(cells[, 1] * cells[, 4] / cells[, 5])
  den <- sum(cells[, 2] * cells[, 3] / cells[, 5])
  num / den
}

rr_mh <- function(tab3d) {
  w <- function(i) {
    a <- tab3d[i, 1, 1]; b <- tab3d[i, 1, 2]
    c <- tab3d[i, 2, 1]; d <- tab3d[i, 2, 2]
    c(a, b, c, d, n = a + b + c + d)
  }
  k <- dim(tab3d)[1]
  cells <- t(vapply(1:k, w, numeric(5)))
  num <- sum(cells[, 1] * (cells[, 3] + cells[, 4]) / cells[, 5])
  den <- sum(cells[, 3] * (cells[, 1] + cells[, 2]) / cells[, 5])
  num / den
}

# same data as tab1, split by sex so both strata have n = 100;
# pooled values below must therefore equal the unstratified tab1 results
tab_sex <- array(0, dim = c(2, 2, 2),
                 dimnames = list(sex = c("male", "female"),
                                 exposed = c("yes", "no"),
                                 disease = c("yes", "no")))
tab_sex["male", , ]   <- matrix(c(15, 35,  8, 42), 2, byrow = TRUE)  # a, b, c, d
tab_sex["female", , ] <- matrix(c(15, 35,  7, 43), 2, byrow = TRUE)  # a, b, c, d
print(tab_sex)

# per stratum: the existing 2x2 functions work unchanged on each slice
round(or2x2(tab_sex["male", , ]), 3)
round(rr2x2(tab_sex["male", , ]), 3)
round(or2x2(tab_sex["female", , ]), 3)
round(rr2x2(tab_sex["female", , ]), 3)

# pooled across strata
round(or_mh(tab_sex), 3)   # MH OR = 2.429, same as tab1 (counts are split evenly)
round(rr_mh(tab_sex), 3)   # MH RR = 2.000

# --- self-check: MH pooling of identical strata equals the unstratified result ---
stopifnot(
  all.equal(unname(round(or_mh(tab_sex), 3)), unname(round(or2x2(tab1)["OR"], 3))),
  all.equal(unname(round(rr_mh(tab_sex), 3)), unname(round(rr2x2(tab1)["RR"], 3)))
)
cat("MH checks passed.\n")

# --- verify against the standard tools on tab1 ---
# glm on long-format data: logit link gives OR, log link gives RR
dat_long <- data.frame(
  disease  = c(rep(1, 30), rep(0, 70), rep(1, 15), rep(0, 85)),
  exposure = c(rep(1, 100), rep(0, 100))
)
exp(coef(glm(disease ~ exposure, binomial, dat_long))[2])   # OR  = 2.429
exp(coef(glm(disease ~ exposure, poisson(link = "log"), dat_long))[2])  # RR = 2.000

# --- self-check: our functions must match the hand-computed values ---
stopifnot(
  all.equal(unname(round(or2x2(tab1)["OR"], 3)), 2.429),
  all.equal(unname(round(rr2x2(tab1)["RR"], 3)), 2.000),
  all.equal(unname(round(or2x2(tab3)["OR"], 3)), 1),
  all.equal(unname(round(rr2x2(tab3)["RR"], 3)), 1)
)
cat("All checks passed.\n")

# Optional: cross-check against epiR (install.packages("epiR"))
# library(epiR)
# epi.2by2(dat = tab1, method = "cohort.count")  # same OR and RR

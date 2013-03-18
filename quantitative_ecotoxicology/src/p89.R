

options(scipen = 1, digits = 5)



## MERCURY <- read.table('p89.csv',
##                       header = TRUE,
##                       sep = ";")



MERCURY <- read.table('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p89.csv', 
                      header = TRUE, 
                      sep = ";")



head(MERCURY)



plot(PCI ~ DAY, data = MERCURY)
abline(v = 22)



# ln-transformation
MERCURY$LPCI <- log(MERCURY$PCI)
# fit linear model for day 31 to 94
mod_slow <- lm(LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
sum_mod_slow <- summary(mod_slow)
sum_mod_slow
exp(coef(mod_slow)[1])



# bias-correction
corr_mod_slow <- exp(sum_mod_slow$sigma^2 / 2) 
# add bias corrected predicts to data.frame
# predict takes the whole data.frame es newdata, so we get predicts for every day.
MERCURY$PRED <- exp(predict(mod_slow, newdata = MERCURY)) * corr_mod_slow
# save C_B and k_B as objects (used later...)
CB <- exp(coef(mod_slow)[1]) * corr_mod_slow
KCB <- abs(coef(mod_slow)[2])



plot(LPCI ~ DAY, data = MERCURY)
points(MERCURY$DAY[1:4], MERCURY$LPCI[1:4], pch = 16)
abline(a = log(CB), b = -KCB)
for(i in 1:4) {
  lines(c(MERCURY$DAY[i], MERCURY$DAY[i]), c(log(MERCURY$PRED[i]), MERCURY$LPCI[i]), lwd = 2)
}



# extract residuals
MERCURY$ERROR  <- MERCURY$PCI - MERCURY$PRED
MERCURY$ERROR[1:4]
# fit linear model to ln(residuals) for Day 3 to 22
mod_fast <- lm(log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY < 22, ])
sum_mod_fast <- summary(mod_fast)
sum_mod_fast
exp(coef(mod_fast)[1])



# bias correction
corr_mod_fast <- exp(sum_mod_fast$sigma^2 / 2) 
# save C_A and k_A as objects
CA <- exp(coef(mod_fast)[1]) * corr_mod_fast
KCA <- abs(coef(mod_fast)[2])



plot(LPCI ~ DAY, data = MERCURY)
abline(mod_slow)
abline(mod_fast, lty = "dotted")
legend("topright", c("slow", "backstripped-fast"), lty=c("solid", "dashed"), cex = 0.8)
# Estimates
c(CA, KCA, CB, KCB)




nls_mod1 <- nls(PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY), 
                data = MERCURY, 
                algorithm = "port",    # to use the bonds
                start = list(KCB = KCB, KCA = KCA, CB = CB, CA = CA),
                lower = c(0, 0, 5000, 20000), 
                upper = c(1, 1, 20000, 45000))
sum_nls_mod1 <- summary(nls_mod1)
sum_nls_mod1



plot(PCI ~ DAY, data = MERCURY, type = "n")
points(MERCURY$DAY, MERCURY$PCI, pch = ifelse(MERCURY$DAY <= 22, 16, 17))
# smooth line
pred_nls_mod1 <- predict(nls_mod1, newdata = data.frame(DAY = seq(0,100, 1)))
lines(seq(0,100, 1), pred_nls_mod1)
legend("topright", c("Fast", "slow"), pch=c(16,17))



sum_mod_slow$sigma^2



exp(sum_mod_slow$sigma^2 / 2)



MERCURY <- read.table('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p89_MERCURY.csv', 
                      header = TRUE, 
                      sep = ";")

MERCURY$LPCI <- log(MERCURY$PCI)
head(MERCURY)

plot(LPCI ~ DAY, data = MERCURY)
abline(v = 22)

mod_slow <- lm(LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
sum_mod_slow <- summary(mod_slow)
# correction
corr_mod_slow <- exp(sum_mod_slow$sigma^2 / 2) 
MERCURY$PRED <- exp(predict(mod_slow, newdata = MERCURY)) * corr_mod_slow
CB <- exp(coef(mod_slow)[1]) * corr_mod_slow
KCB <- abs(coef(mod_slow)[2])



MERCURY$ERROR  <- MERCURY$PCI - MERCURY$PRED
mod_fast <- lm(log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY <= 22, ])
# na.omit for negative value
sum_mod_fast <- summary(mod_fast)
corr_mod_fast <- exp(sum_mod_fast$sigma^2 / 2) 
CA <- exp(coef(mod_fast)[1]) * corr_mod_fast
KCA <- abs(coef(mod_fast)[2])

plot(LPCI ~ DAY, data = MERCURY)
abline(mod_slow)
abline(mod_fast, lty = "dotted")
legend("topright", c("slow", "backstripped-fast"), lty=c("solid", "dashed"), cex = 0.8)



nls_mod1 <- nls(PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY), 
                data = MERCURY, 
                algorithm = "port",    # to use the bonds
                start = list(KCB = KCB, KCA = KCA, CB = CB, CA = CA),
                lower = c(0, 0, 5000, 20000), 
                upper = c(1, 1, 20000, 45000))

plot(PCI ~ DAY, data = MERCURY, type = "n")
points(MERCURY$DAY, MERCURY$PCI, pch = ifelse(MERCURY$DAY <= 22, 16, 17))
# smooth line
pred_nls_mod1 <- predict(nls_mod1, newdata = data.frame(DAY = seq(0,100, 1)))
lines(seq(0,100, 1), pred_nls_mod1)
legend("topright", c("Fast", "slow"), pch=c(16,17))

summary(nls_mod1)
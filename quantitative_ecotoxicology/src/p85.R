OYSTERZN <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p85.csv", 
                      header = TRUE, 
                      sep = ";")

plot(ZINC ~ DAY, data = OYSTERZN)

# without weights
mod <- nls(ZINC ~ INITACT * exp((-(KE+0.00283))*DAY), OYSTERZN,
           start = list(KE=0.004, INITACT = 500))
summary(mod)
pred_mod <- predict(mod)
res_mod <- resid(mod)
# Residuals
plot(OYSTERZN$DAY, res_mod)
abline(h = 0)


# with fixed weights
mod2 <- nls(ZINC ~ INITACT * exp((-(KE+0.00283))*DAY), OYSTERZN,
            weights = OYSTERZN$WNLIN, 
            start = list(KE=0.004, INITACT = 500))
summary(mod2)
pred_mod2 <- predict(mod2)

# Data + fitted
plot(ZINC ~ DAY, data = OYSTERZN)
lines(OYSTERZN$DAY, pred_mod)
lines(OYSTERZN$DAY, pred_mod2, lty = "dashed")



# Lm
OYSTERZN$LZINC <- log(OYSTERZN$ZINC)
mod3 <- lm(LZINC ~ DAY, data = OYSTERZN)
summary(mod3)
pred_mod3 <- predict(mod3)

# data + fitted
plot(LZINC ~ DAY, OYSTERZN)
lines(OYSTERZN$DAY, pred_mod3)

# residuals
plot(mod3, which = 1)

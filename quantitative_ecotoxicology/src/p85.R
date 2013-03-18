

OYSTERZN <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p85.csv", 
                      header = TRUE, 
                      sep = ";")




## OYSTERZN <- read.table("p85.csv",
##                       header = TRUE,
##                       sep = ";")



head(OYSTERZN)



plot(ZINC ~ DAY, data = OYSTERZN)



mod_noweight <- nls(ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY), 
           data = OYSTERZN,
           start = list(KE = 0.004, INITACT = 500))



summary(mod_noweight)



sum_noweight <- summary(mod_noweight)



res_mod_noweight <- resid(mod_noweight)
plot(OYSTERZN$DAY, res_mod_noweight)
abline(h = 0)



OYSTERZN$WNLIN <- OYSTERZN$DAY^2



mod_weight <- nls(ZINC ~ INITACT * exp((-(KE+0.00283))*DAY), 
                  data = OYSTERZN,
                  weights = OYSTERZN$WNLIN, 
                  start = list(KE = 0.004, INITACT = 500))
summary(mod_weight)



sum_weight <- summary(mod_weight)



# extract the fitted values
fit_mod_noweight <- fitted(mod_noweight)
fit_mod_weight <- fitted(mod_weight)
# plot data
plot(ZINC ~ DAY, data = OYSTERZN)
# add fitted values
lines(OYSTERZN$DAY, fit_mod_noweight)
lines(OYSTERZN$DAY, fit_mod_weight, lty = "dashed")
# add legend
legend("topright", legend = c("nonweighted", "weighted"), lty=c("solid", "dashed"))



OYSTERZN$LZINC <- log(OYSTERZN$ZINC)



plot(LZINC ~ DAY, OYSTERZN)



mod_lm <- lm(LZINC ~ DAY, data = OYSTERZN)
sum_lm <- summary(mod_lm)
sum_lm



# fitted values
fit_mod_lm <- fitted(mod_lm)

# data + fitted
plot(LZINC ~ DAY, OYSTERZN)
lines(OYSTERZN$DAY, fit_mod_lm)

# residuals
plot(mod_lm, which = 1)



# MSE
sum_lm$sigma^2



# unbiased estimate for C_0: exp(MSE/2) * exp(Intercept)
exp(sum_lm$sigma^2 / 2) * exp(sum_lm$coefficients[1, 1])



## sum_lm$coefficients[1, 1]



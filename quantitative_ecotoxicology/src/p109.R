
options(scipen = 1, digits = 5)



require(knitr)
opts_chunk$set(out.width="400px", fig.height=6, fig.width=6)



require(RCurl)
# Accumulation
url_accum <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_accum.csv",
ssl.verifypeer = FALSE)
ACCUM <- read.table(text = url_accum, header = TRUE, sep = ";")




head(ACCUM)



mod_accum <- nls(BRPHOS ~ KU / KE * 10.5 * (1 - exp(-KE * HOUR)),
           data = ACCUM, 
           start = list(KU = 100, KE = 0.01))



summary(mod_accum)



HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1) 
# Raw data
plot(ACCUM, main = 'Accumulation')
# add model
lines(HOUR_pred, predict(mod_accum, newdata = data.frame(HOUR = HOUR_pred)))



# Elimination data
url_elimin <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_elimin.csv")
ELIMIN <- read.table(text = url_elimin, header = TRUE, sep = ";")



head(ELIMIN)
plot(ELIMIN)



ELIMIN$LBROMO <- log(ELIMIN$BRPHOS)



mod_elimin_lm <- lm(LBROMO ~ HOUR, data = ELIMIN)
summary(mod_elimin_lm)



par(mfrow = c(1, 2))
# plot linearized model
plot(LBROMO ~ HOUR, data = ELIMIN, main = 'Data + Model')
# add regression line
abline(mod_elimin_lm)
# plot residuals
plot(fitted(mod_elimin_lm), residuals(mod_elimin_lm), main = 'Residuals')
abline(h = 0, lty = 'dotted')



mod_accum2 <- nls(BRPHOS ~ KU / -coef(mod_elimin_lm)[2] * 10.5 * (1 - exp(coef(mod_elimin_lm)[2] * HOUR)),
           data = ACCUM, 
           start = list(KU = 100))
summary(mod_accum2)



par(mfrow=c(1,2))
HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1) 
# Raw data
plot(ACCUM, main = 'Accumulation')
# add model
lines(HOUR_pred, predict(mod_accum2, newdata = data.frame(HOUR = HOUR_pred)))
plot(fitted(mod_accum2), residuals(mod_accum2))



# Residual sum of squares
mod_accum$m$deviance()
mod_accum2$m$deviance()
# AIC
AIC(mod_accum)
AIC(mod_accum2)



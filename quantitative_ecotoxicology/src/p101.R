
options(scipen = 1, digits = 5)



require(knitr)
opts_chunk$set(out.width="400px", fig.height=6, fig.width=6)



require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p101.csv",
ssl.verifypeer = FALSE)
ZINC <- read.table(text = url, header = TRUE, sep = ";")



head(ZINC)



mod_nls <- nls(N ~ (K*C*M)/(1+K*C), data = ZINC, 
           start = list(K = 3, M = 9), lower = 0, 
           algorithm = 'port')



summary(mod_nls)



plot(ZINC$C, ZINC$N, xlab = 'C', ylab = 'N')
# generate C-values to predict
x_n <- seq(min(ZINC$C), max(ZINC$C), length.out=200)
# add predicts to plot
lines(x_n, predict(mod_nls, newdata = data.frame(C = x_n)))



ZINC$Y <- ZINC$C / ZINC$N



mod_lm <- lm(Y ~ C, data = ZINC)
plot(ZINC$C, ZINC$Y, ylab = 'C/N', xlab = 'C')
abline(mod_lm)
summary(mod_lm)



ZINC$WGT = ZINC$N^4 / ZINC$C^2



mod_wgt <- lm(Y ~ C, data = ZINC, weights = ZINC$WGT)
summary(mod_wgt)



coef(mod_wgt)[2] / coef(mod_wgt)[1]



1 / coef(mod_wgt)[2]



par(mfrow = c(1,2))
# lm
plot(mod_lm, which = 1, main='linear model without weights')
# nls
plot(fitted(mod_nls), residuals(mod_nls), xlab = 'fitted', ylab = 'Residuals', main = 'nonlinear regression')
abline(h = 0, lty = 'dotted')



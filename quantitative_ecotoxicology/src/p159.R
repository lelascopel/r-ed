
options(scipen = 1, digits = 5)
require(knitcitations)
cite_options(linked=TRUE)



require(knitr)
opts_chunk$set(fig.height=6, fig.width=6)



require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p160.csv",
ssl.verifypeer = FALSE)
NAP <- read.table(text = url, header = TRUE, sep = ";")



head(NAP)



NAP$PROP <- NAP$DEAD / NAP$TOTAL



plot(NAP$CONC, NAP$PROP, 
     pch = 16, 
     xlab = expression(paste('Concentration (', mu, 'g/L)')),
     ylab = 'Proportion Dead',
     main = 'Raw data')



contr_m <- t.test(NAP$PROP[NAP$CONC==0])
contr_m



## extract the values from t.test-object
# mean
contr_m$estimate
# CI
contr_m$conf.int



d_control <- mean(NAP$PROP[NAP$CONC == 0])
d_control



NAP$PROP_c <- (NAP$PROP - d_control) / (1 - d_control)
NAP$PROP_c



require(drc)
mod1 <- drm(PROP ~ CONC, data = NAP, fct = LL.2())



mselect(mod1, fctList= list(LL.3(), LL.4(), LL.5(), W1.2(), W1.3(), W1.4()))



plot(mod1, broken = TRUE, type = 'all', bp = 500, xt = seq(500,3000,500))
mtext('Dose-Response-Model - LL2.2', 3)



NAP$PROP_c[NAP$PROP_c < 0 | NAP$CONC == 0] <- 0



mod2 <- drm(PROP_c ~ CONC, data = NAP, fct = LL.2())



mselect(mod2, fctList= list(LL.3(), LL.4(), LL.5(), W1.2(), W1.3(), W1.4()))
mod2 <- update(mod2, fct = W1.2())



plot(mod2, broken = TRUE, type = 'all', bp = 500, xt = seq(500,3000,500))
mtext('Corrected mortalities - W1.2', 3)



mod3 <- drm(PROP ~ CONC, data = NAP, fct = LL.3u())
plot(mod3, broken = TRUE, type = 'all', bp = 500, xt = seq(500,3000,500))
mtext('Free (estimated) lower limit - LL3.u', 3)



summary(mod3)



mselect(mod3, fctList = list(LL.2(), LL2.3u()))



ED(mod1, 50, interval='delta')
ED(mod2, 50, interval='delta')
ED(mod3, 50, interval='delta')



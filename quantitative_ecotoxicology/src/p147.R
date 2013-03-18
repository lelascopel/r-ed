
options(scipen = 1, digits = 5)
require(knitcitations)
cite_options(linked=TRUE)



require(knitr)
opts_chunk$set(out.width="400px", fig.height=6, fig.width=6)



require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p147.csv",
ssl.verifypeer = FALSE)
salt <- read.table(text = url, header = TRUE, sep = ";")



head(salt)



salt$prop <- salt$DEAD / salt$TOTAL



plot(salt$CONC, salt$prop, xlab = 'Concentration', ylab = 'Proportion dead', log='x')



require(drc)
mod <- drm(prop ~ CONC, data = salt, fct =  LL.2())



mselect(mod, fctList=list(W1.2(),G.2()))



# raw data
plot(prop ~ CONC, data = salt, xlim = c(9,21), ylim = c(0,1), ylab = 'Proportion dead', xlab = 'Concentration')

conc_pred <- seq(9, 21, 0.1)
lines(conc_pred, predict(mod, newdata = data.frame(CONC = conc_pred)))



ED(mod, 50, interval='delta')



bibliography('html')



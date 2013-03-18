
options(scipen = 1, digits = 5)



require(knitr)
opts_chunk$set(out.width="400px", fig.height=6, fig.width=6)



require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p108.csv",
ssl.verifypeer = FALSE)
MERCURY <- read.table(text = url, header = TRUE, sep = ";")



head(MERCURY)



plot(MERCURY)



mod <- nls(HG ~ KU / KE * 0.24 * (1 - exp(-KE * DAY)), 
           data = MERCURY, 
           start = list(KU = 1000, KE = 0.5))



summary(mod)



BCF = coef(mod)[1] / coef(mod)[2]
BCF



BCF * 0.24



DAY_pred <- seq(0, 6, by = 0.1) 
# Raw data
plot(MERCURY)
# add model
lines(DAY_pred, predict(mod, newdata = data.frame(DAY = DAY_pred)))
# add model-equation
text(3, 100, bquote(HG == .(BCF*0.24)%.%(1-exp(-.(coef(mod)[2])%.%DAY))))




options(scipen = 1, digits = 5)
require(knitcitations)
cite_options(linked=TRUE)



require(knitr)
opts_chunk$set(fig.height=6, fig.width=6)



require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/TOXICITY.csv",
ssl.verifypeer = FALSE)
TOXICITY <- read.table(text = url, header = TRUE)
head(TOXICITY)
summary(TOXICITY)



TOXICITY$FLAG <- ifelse(TOXICITY$TTD > 96, 1, 2)



require(survival)
mod <- survfit(Surv(TTD, FLAG) ~ PPT + strata(TANK), data = TOXICITY)
plot(mod, col = rep(1:7, each=2), mark.time=FALSE)
legend('bottomleft', legend = sort(unique(TOXICITY$PPT)), col=1:7, lty = 1)



survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==10.3, ], rho = 0)
survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==10.8, ], rho = 0)
survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==11.6, ], rho = 0)
survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==13.2, ], rho = 0)
survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==15.8, ], rho = 0)



for(i in sort(unique(TOXICITY$PPT)[-c(2,7)])) {
  cat('\n', i, '\n')
  print(survdiff(Surv(TTD, FLAG) ~ TANK, data = TOXICITY[TOXICITY$PPT==i, ], rho = 1))
}



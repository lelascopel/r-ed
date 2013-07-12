
options(scipen = 1, digits = 5)
require(knitcitations)
cite_options(linked=TRUE)



require(knitr)
opts_chunk$set(fig.height=6, fig.width=6)



# download data from github
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/TOXICITY.csv", ssl.verifypeer = FALSE)
TOXICITY <- read.table(text = url, header = TRUE)
# add status
TOXICITY$FLAG <- ifelse(TOXICITY$TTD > 96, 1, 2)



require(survival)
# fit different models
mod_exp <- survreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'exponential')
mod_wei <- survreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'weibull')
mod_lnorm <- survreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'lognorm')
mod_llog <- survreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'loglogistic')



# extract loglik and save into data.frame
df <- data.frame(model = c('exp', 'weibull', 'lnorm', 'loglog'),
                 logLik = c(mod_exp$loglik[2],
                            mod_wei$loglik[2],
                            mod_lnorm$loglik[2],
                            mod_llog$loglik[2]))



# add AIC as column to data.frame
df$AIC <- c(extractAIC(mod_exp)[2],
            extractAIC(mod_wei)[2],
            extractAIC(mod_lnorm)[2],
            extractAIC(mod_llog)[2])
df



summary(mod_wei)



# newdata, all combinations of WETWT and PPT
newtox <- expand.grid(WETWT = seq(0, 1.5, length.out=100),
            PPT = seq(12.5, 20, by = 2.5))
# predict ttd for newdata
newtox$preds <- predict(mod_wei, newdata=newtox, type='quantile', p = 0.5)

# plot
require(ggplot2)
ggplot(newtox, aes(x = WETWT, y = preds, col = factor(PPT), group = PPT)) +
  geom_line() +
  theme_bw() +
  labs(x = 'Wet Weight (g)', y = 'Median-Time-To-Death (h)', col = 'Concentration (g/L)') +
  coord_cartesian(ylim=c(0,100))



# plot KM
km <- survfit(Surv(TTD, FLAG) ~ PPT, data = TOXICITY)
plot(km, col = 1:7)

# Fit model only for Concentrations
mod_ppt <- survreg(Surv(TTD, FLAG) ~ PPT, data = TOXICITY, dist = 'weibull')

# add model to plot
PPT_u <- sort(unique(TOXICITY$PPT))
for(i in seq_along(PPT_u)){
  lines(predict(mod_ppt, 
                newdata = list(PPT=PPT_u[i]),
                type = "quantile",
                p = 1:99/100),
                99:1/100, 
                col=i)
}



# opposite of %in% : 'not in'
`%!in%` <- Negate(`%in%`)

# remove 0 and 20.1 treatments
TOXICITY_sub <- TOXICITY[TOXICITY$PPT %!in% c(0, 20.1), , drop = TRUE]
# kaplan-meier estimates
km <- survfit(Surv(TTD, FLAG) ~ PPT, data = TOXICITY_sub)
fac <- factor(rep(sort(unique(TOXICITY_sub$PPT)), km$strata))
cols <- 1:5
# plot tranformed time vs. survival
plot(log(km$time), log(-log(km$surv)), col=cols[fac], pch = 16)
# add legend
legend('bottomright', legend=levels(fac), col=cols, pch = 16, cex = 0.7)



#### Define loglogistic distribution for flexsurvreg
library(eha)  ## make "dllogis" and "pllogis" available to the working environment 
custom.llogis <- list(name="llogis",
                      pars=c("shape","scale"),
                      location="scale",
                      transforms=c(log, log),
                      inv.transforms=c(exp, exp),
                      inits=function(t){ c(1, median(t)) })



require(flexsurv)
# fit all models using flexsurvreg
mod_flex_exp <- flexsurvreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'exp')
mod_flex_wei <- flexsurvreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'weibull')
mod_flex_lnorm <- flexsurvreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'lnorm')
mod_flex_llog <- flexsurvreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = custom.llogis)
mod_flex_ggamma <- flexsurvreg(Surv(TTD, FLAG) ~ PPT + WETWT, data = TOXICITY, dist = 'gengamma')



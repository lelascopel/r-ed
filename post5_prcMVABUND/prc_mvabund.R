setwd("/home/edisz/Documents/Uni/Projects/blog/post5_prcMVABUND/")

require(vegan)
require(mvabund)


data(pyrifos)
# The data has been log-transformed, we want to use the raw abundances:
pyrifos <- round((exp(pyrifos) - 1)/10)

week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))

# mvabund is very computer intensiv, so I 
# take only the 20 most abundant species in week -1
# to reduce the computation time for destration
take <- pyrifos[ , order(apply(pyrifos[week == -1, ], 2, max), decreasing=TRUE)[1:10]]

abudat <- mvabund(take)

# Global test
mod1 <- manyglm(abudat ~ week * dose, family="n")
plot(mod1)
tt <- system.time(aov_mod1 <- anova(mod1))
aov_mod1
saveRDS(aov_mod1, "aov_mod1.rds")
saveRDS(tt, "tt.rds")

# test per week, would be nicer with contrasts
out <- NULL
for(i in levels(week)) {
  take_spec <- abudat[week == i, ]
  take_dose <- dose[week == i]
  mod <- manyglm(take_spec ~ take_dose)
  out[[i]]$summary <- summary(mod, p.uni = "adjusted")
  out[[i]]$fitted <- fitted(mod)
}
out
saveRDS(out, "out.rds")



###
## plot 
require(ggplot2)
p <- ggplot(df, aes(x = .id, y = lr, col = dose)) +
  geom_line() +
  ylab("Likelihood-Ratio") +
  xlab("Week")
p


## anderer Plot
require(reshape2)
ftd <- predict(mod1, type = "response")   # predicts
# bring to long format
pyr_backtransformed_melt <- melt(cbind(week, dose, pyr_backtransformed))
ftd_melt <- melt(data.frame(week = week, dose = dose, ftd))

# merge fitted and raw
merged <- merge(pyr_backtransformed_melt, ftd_melt, by = c("week", "dose", "variable"))
merged$week <- as.numeric(as.character(merged$week))
# Eventuell noch confidence-Intervalle rein...

ggplot(merged) +
  geom_point(aes(x = week, y = value.x, col = dose)) +
  facet_wrap(~variable, scales = "free_y") +
  geom_line(aes(x = week, y = value.y, col = dose)) +
  #  geom_smooth(aes(x = week, y = value.y, col = dose, ymin = lcl, ymax = ucl), stat="identity") +
  #  scale_y_log10() +
  geom_vline(aes(xintercept = 0), col ="red") +
  theme_bw()

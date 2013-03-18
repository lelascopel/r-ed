
options(scipen = 1, digits = 5)



## LEAD <- read.table('p94.csv',
##                       header = TRUE,
##                       sep = ";")



LEAD <- read.table('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p94.csv', 
                      header = TRUE, 
                      sep = ";")



head(LEAD)



plot(LEAD ~ DAY, LEAD)



LEAD$LLEAD <- log(LEAD$LEAD)
LEAD$LDAY <- log(LEAD$DAY)
plot(LLEAD ~ LDAY, LEAD)



# fit model
mod <- lm(LLEAD ~ LDAY, data = LEAD)



plot(mod, which = 1)



mod_sum <- summary(mod)
mod_sum



plot(LLEAD ~ LDAY, LEAD)
abline(mod)



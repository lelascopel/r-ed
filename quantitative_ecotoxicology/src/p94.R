LEAD <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S94.csv"), 
                      header = TRUE, 
                      sep = ";")

LEAD$LLEAD <- log(LEAD$LEAD)
LEAD$LDAY <- log(LEAD$DAY)

mod <- lm(LLEAD ~ LDAY, data = LEAD)
plot(LLEAD ~ LDAY, data = LEAD)
abline(mod)
summary(mod)
# diagnostics
par(mfrow = c(2,2))
plot(mod)
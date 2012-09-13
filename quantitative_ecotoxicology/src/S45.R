CADMIUM <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S45.csv"), 
                      header = TRUE, 
                      sep = ";")


CADMIUM$FLIP <- abs(CADMIUM$CD - 100)

require(survival)
# log-rank
fit <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
fit
# Peto & Peto modification of the Gehan-Wilcoxon test
fit2 <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
fit2
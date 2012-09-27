CADMIUM <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p45.csv", 
                  header = TRUE, 
                  sep = ";")


CADMIUM$FLIP <- abs(CADMIUM$CD - 100)
CADMIUM <- CADMIUM[order(CADMIUM$SITE, CADMIUM$FLIP), ]

require(survival)
# log-rank
fit <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
fit
# Peto & Peto modification of the Gehan-Wilcoxon test
fit2 <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
fit2
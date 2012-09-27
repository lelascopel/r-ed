SULFATE <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p39.csv", 
                  header = TRUE, 
                  sep = ";")

SULFATE$FLIP <- abs(SULFATE$SO4 - 8)

require(survival)
fit <- survfit(Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type="plain") 
fit
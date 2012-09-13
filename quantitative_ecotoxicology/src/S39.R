SULFATE <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/SULFATE.csv"), 
                  header = TRUE, 
                  sep = ";")

SULFATE$FLIP <- abs(SULFATE$SO4 - 8)

require(survival)
fit <- survfit(Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type="plain") 
fit
SULFATE <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S42.csv"), 
                      header = TRUE, 
                      sep = ";")

wilcox.test(SO4 ~ SITE, SULFATE, correct = TRUE, alternative = "greater")
wilcox.test(SO4 ~ SITE, SULFATE, correct = TRUE)
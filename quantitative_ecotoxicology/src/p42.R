

SULFATE <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p42.csv", 
                  header = TRUE, 
                  sep = ";")



## SULFATE <- read.table("p42.csv",
##                   header = TRUE,
##                   sep = ";")



require(ggplot2)
ggplot(SULFATE, aes(x = SITE, y = SO4)) +
  geom_violin()



wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE)



wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE, alternative = "greater")



require(coin)
wilcox_test(SO4 ~ SITE, SULFATE, distribution="exact", conf.int = TRUE)



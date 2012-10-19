QCBLKS<- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p47_QCBLKS.csv", 
                      header = TRUE, 
                      sep = ";")

LEADBLK <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p47_LEADBLK.csv", 
                      header = TRUE, 
                      sep = ";")

require(qcc)
a <- qcc(df_c, "xbar")
plot(a)
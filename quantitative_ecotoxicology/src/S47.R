QCBLKS<- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S47_QCBLKS.csv"), 
                      header = TRUE, 
                      sep = ";")

LEADBLK <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S47_LEADBLK.csv"), 
                      header = TRUE, 
                      sep = ";")

require(qcc)
a <- qcc(df_c, "xbar")
plot(a)


SO4 <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p35.csv", 
                  header = TRUE, 
                  sep = ";")



## SO4 <- read.table("p35.csv",
##                  header = TRUE,
##                  sep = ";")



SO4$rem2 <- ifelse(SO4$rem == "<", TRUE, FALSE)
SO4



require(NADA)
rs <- ros(SO4$value, SO4$rem2)
print(rs)



SO4 <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/SO4.csv"), 
                  header = TRUE, 
                  sep = ";")
SO4$rem <- ifelse(SO4$rem == "<", TRUE, FALSE)

require(NADA)
ros(SO4$value, SO4$rem)
ALL <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/ALL.csv"), 
                 header = TRUE, 
                 sep = ";")

winsori <- function (x, width = 2)
{
  if(is.unsorted(x))
    stop("Values must be sorted!")
  n <- length(x)
  x[1:width] <- x[width + 1]
  x[(n-width):n] <- x[(n-width)]
  x
}

ALL$SO4_win <- winsori(ALL$SO4)
mean(ALL$SO4_win)
sd(ALL$SO4_win)
sw <- function(x, width = 2){
  n <- length(x)
  sd(x) * (n - 1) / (n - 2*width -1)
}
sw(ALL$SO4_win)
mean(ALL$SO4, trim=2/21)


ALL <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S33.csv"), 
                 header = TRUE, 
                 sep = ";")
head(ALL)

winsori <- function (x, width = 2)
{
  # check if sorted
  if(is.unsorted(x))
    stop("Values must be sorted!")
  # get number of observations
  n <- length(x)
  # Replace lowest
  x[1:width] <- x[width + 1]
  # replace highest
  x[(n - width + 1):n] <- x[(n-width)]
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
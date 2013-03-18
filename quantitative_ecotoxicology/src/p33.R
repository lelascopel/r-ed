

ALL <- read.table("/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/p33.csv", 
                 header = TRUE, 
                 sep = ";")



## ALL <- read.table("p33.csv",
##                  header = TRUE,
##                  sep = ";")



str(ALL)
ALL$SO4



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
# display the first and 5 last rows
ALL[c(1:5, 17:21), ]



# mean
mean(ALL$SO4_win)
# standard deviation
sd(ALL$SO4_win)



sw <- function(x, width = 2){
  n <- length(x)
  sd(x) * (n - 1) / (n - 2*width -1)
}
sw(ALL$SO4_win)



mean(ALL$SO4, trim=2/21)



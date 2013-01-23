## function to read data from SAS/pdf
## start copying from DATALINES till ';'

read_new <- function(x){
  con <- file(x)
  tmp <- readLines(con) # Read one line 
  close(con)
  
  require(stringr)
  var <- strsplit(tmp[1], split = " ")
  var_a <- which(var[[1]] == "@@;")
  vars <- var[[1]][2:(var_a - 1)]
  nvar <- length(vars)
  
  data_start <- grep("DATALINES;", tmp)
  data <- tmp[(data_start + 1):(length(tmp)-1)]
  data <- paste(data, collapse= " ")
  data <- unlist(strsplit(data, split=" "))
  
  out <- list()
  for(i in seq_along(vars)){
    se <- 1:nvar
    take <- rep(se, times = length(data) / nvar)    
    out[[vars[i]]] <- data[take == i]
  }
  out <- data.frame(out)
  out
}


files <- list.files('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/data/', full.names=TRUE, pattern="*.txt")
for(i in seq_along(files)){
  print(files[i])
  df <- read_new(files[i])
  write.table(df, gsub(".txt", ".csv", files[i]), row.names=FALSE, sep = ";")
}



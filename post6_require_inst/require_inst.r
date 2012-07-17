require_inst <- function(x) {
  if(!require(x, character.only = TRUE)) {
    cat("\nPackage not installed!\nEnter 'y' to install package, 'n' to exit:\n") 
    y <- readline()
    if(y != "y") {
      stop("Aborted.")
    } else {
      install.packages(x)
      if(!require(x, character.only = TRUE)){
        stop("Could not install package!")
      }
    }
  } 
  else {
    require(x, character.only = TRUE)
  }
}

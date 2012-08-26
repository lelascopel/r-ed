require(RCurl)
require(XML)

get_uid <- function(x){
  x <- gsub(" ", "+", x)
  searchurl <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=",
                     x, sep = "")
  tt <- getURL(searchurl)
  ttp <- xmlTreeParse(tt, useInternalNodes = TRUE) 
  res <- xpathSApply(ttp, "//eSearchResult/IdList/Id", xmlValue)
  # if xpath is not found return NA
  if(length(res) == 0) { 
    out <- NA
  } 
  else {
    out <- res
  }
  #NCBI limits requests to three per second
  Sys.sleep(0.33)
  return(out)
}

uid <- get_uid("Hydropsyche angustipennis")
uid

get_classification <- function(x){
  baseurl <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"
  ID <- paste("ID=", x, sep ="")
  searchurl <- paste(baseurl, ID, sep = "&")
  tt <- getURL(searchurl)
  ttp <- xmlTreeParse(tt, useInternalNodes = TRUE)
  out <- data.frame(ScientificName = xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/ScientificName", xmlValue),
                    Rank = xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/Rank", xmlValue),
                    UID = xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/TaxId", xmlValue))
  #NCBI limits requests to three per second
  Sys.sleep(0.33)
  return(out)
}
get_classification(uid)
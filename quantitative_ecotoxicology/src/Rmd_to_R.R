files <- list.files('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/rmd', pattern='.Rmd')
infiles <- file.path('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/rmd', files)
out <- gsub('.Rmd', '.R', files)
outfiles <- file.path('/home/edisz/Documents/Uni/Projects/blog/quantitative_ecotoxicology/src', out)

for(i in seq_along(infiles)){
  purl(input = infiles[i], output = outfiles[i], documentation = 0)
}

options(scipen = 1, digits = 5)
require(knitcitations)
cite_options(linked=TRUE)



require(knitr)
opts_chunk$set(fig.height=6, fig.width=6)



TEST <- matrix(c(1,19,6,14), byrow=TRUE, ncol = 2, 
               dimnames=list(c('Tank_A', 'Tank_B'), c('Number_Dead', 'Number_Surviving')))
TEST



fisher.test(TEST)
fisher.test(TEST, alternative='greater')
fisher.test(TEST, alternative='less')



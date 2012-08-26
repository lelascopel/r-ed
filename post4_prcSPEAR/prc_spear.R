setwd("/home/edisz/Documents/Uni/Projects/blog/post4_prcSPEAR/")

require(vegan)    # for the data
require(xlsx)     # writing to .xls
require(RCurl)    # to read .csv from github
require(reshape2) # tranform data formats (long / wide)
require(gdata)    # read.xls (don't know why read.xlsx doesn't work)
require(ggplot2)  # graphics
require(plyr)     # aggregate and transform data

# pyrifos data
data(pyrifos)
# The data has been log-transformed, we want to use the raw abundances:
pyrifos <- round((exp(pyrifos) - 1)/10)

pyrifos$week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
pyrifos$dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
pyrifos$ditch <- gl(12, 1, length=132)


# bring data to long format
pyrifos_melt <- melt(pyrifos, id = c("dose", "week", "ditch"))

# get lookup table for taxa-names from my github repository
csv <- getURL("https://raw.github.com/EDiLD/r-ed/master/post4_prcSPEAR/taxa_names.csv")
taxa_names <- read.table(textConnection(csv), 
                         header = TRUE, sep = ";")

# replace taxa names
pyrifos_melt$variable <- taxa_names$taxa[match(pyrifos_melt$variable, taxa_names$abbrv)]

# remove unmatched names
pyrifos_melt <- pyrifos_melt[!is.na(pyrifos_melt$variable), ]

# write to xls
write.xlsx2(pyrifos_melt, 
            file = "spear_in.xlsx", 
            row.names = FALSE)

#### 
####
####

# read xls from SPEAR-calculator
spear <- read.xls("spear_out.xls", 
                  sheet = 3)
spear$Treatment <- factor(spear$Treatment)

# plot raw
p <- ggplot(spear, aes(x = Time.Point, y = SPEARpesticides, col = Treatment)) +
  geom_point() +
  geom_line(aes(group = Replicates)) +
  xlab("Weeks")
p

# plot means
spear_means <- ddply(spear, .(Treatment, Time.Point), summarise, 
                     mean_spear = mean(SPEARpesticides))

p2 <- ggplot(spear_means, aes(x = Time.Point, y = mean_spear, col = Treatment)) +
  geom_point() +
  geom_line() +
  xlab("Weeks") +
  ylab("SPEARpesticides")
p2

# plot difference to control
spear_relcont <- ddply(spear_means, .(Time.Point), transform, 
      mean_spear = mean_spear - mean_spear[Treatment == 0])

p3 <- ggplot(spear_relcont, aes(x = Time.Point, y = mean_spear, col = Treatment)) +
  geom_point() +
  geom_line() +
  xlab("Weeks") +
  ylab(expression(paste(Delta, "SPEARpesticides")))
p3


##### Check ylab! with beketov.





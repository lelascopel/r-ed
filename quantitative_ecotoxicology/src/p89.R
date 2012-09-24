MERCURY <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S89_MERCURY.csv"), 
                      header = TRUE, 
                      sep = ";")

MERCURY$LPCI <- log(MERCURY$PCI)
EXPO <- lm(LPCI ~ DAY, data = MERCURY)

PEXPO <- fitted(EXPO)
plot(LPCI ~ DAY, data = MERCURY)
abline(EXPO)
plot(EXPO, which = 1)

CB <- exp(coef(EXPO)[1])
KCB <- -coef(EXPO)[2]

##############################################
FAST <- read.table(file.path(getwd(), "quantitative_ecotoxicology/data/S89_FAST.csv"), 
                      header = TRUE, 
                      sep = ";")

plot(FAST)
names(FAST)[2] <- "PCI"
FAST$LPCI <- log(FAST$PCI)  # no log of neg number!
FEXPO <- lm(LPCI ~ DAY, data = FAST)

PFEXPO <- fitted(FEXPO)
plot(LPCI ~ DAY, data = FAST)
abline(FEXPO)

plot(FEXPO, which = 1)

CA <- exp(coef(FEXPO)[1])
KCA <- -coef(FEXPO)[2]

######################################

nls_mod1 <- nls(PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY), 
                data = MERCURY, 
                algorithm = "port",    # to use the bonds
                start = list(KCB = KCB, KCA = KCA, CB = CB, CA = CA),
                lower = c(0, 0, 5000, 20000), 
                upper = c(1, 1, 20000, 45000))

plot(PCI ~ DAY, data = MERCURY)
lines(MERCURY$DAY, fitted(nls_mod1))

summary(nls_mod1)

### Quantitative Ecotoxicology with R

There is a fresh book about [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647):
> Newman, Michael C. (2012). Quantitative Ecotoxicology, Second Edition. CRC Press


Unfortunately all examples are written in SAS, to which I have no access (5.530€ for a licence of SAS Analytics Pro).

Since I am a poor student who likes R and ecotoxicology, I will try to reproduce the examples given in the book.


All code and data are available at my [github-repository](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology).
### Quantitative Ecotoxicology, Page 101, Example 3.6, Langmuir

This is example 3.6 on page 101 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about adsorption and how to fit an adsorption model to data.

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p101.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p101.csv",
ssl.verifypeer = FALSE)
ZINC <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(ZINC)
```

```
##      N     C
## 1 0.75 0.030
## 2 1.40 0.069
## 3 1.95 0.118
## 4 2.51 0.166
## 5 3.03 0.217
## 6 3.53 0.270
```


So we have a data.frame with two columns,
where N = amount adsorbed (mmol) per unit mass (g) and  C = equilibrium concentration in the aqueous phase (mmol/ml).

We want fit a Langmuir Model (Equation 3.28 in the book) to this data. 

The three methods described are:

* Nonlinear Regression
* linear transformation
* linear transformation with weighting



#### Nonlinear Regression

```r
mod_nls <- nls(N ~ (K * C * M)/(1 + K * C), data = ZINC, start = list(K = 3, 
    M = 9), lower = 0, algorithm = "port")
```

This fits the model 

$$ N = \frac{KCM}{1+KC} $$ 

to the data. 

We supplied some starting values and specified the lower bonds for K and M as 0 (bonds can only be used with the port algorithm).

This gives us the estimates for K and M as:

```r
summary(mod_nls)
```

```
## 
## Formula: N ~ (K * C * M)/(1 + K * C)
## 
## Parameters:
##   Estimate Std. Error t value Pr(>|t|)    
## K    2.097      0.188    11.1  3.8e-06 ***
## M    9.899      0.521    19.0  6.1e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0929 on 8 degrees of freedom
## 
## Algorithm "port", convergence message: relative convergence (4)
```


* $K = 2.097 \pm 0.188$
* $M = 9.899 \pm 0.521$

The t and p-values of this output are not of interest for us (tests if the parameters deviate from 0).

We can plot the raw data and the model easily using the predict-function:

```r
plot(ZINC$C, ZINC$N, xlab = "C", ylab = "N")
# generate C-values to predict
x_n <- seq(min(ZINC$C), max(ZINC$C), length.out = 200)
# add predicts to plot
lines(x_n, predict(mod_nls, newdata = data.frame(C = x_n)))
```

<img src="figure/plot-nls.png" title="plot of chunk plot-nls" alt="plot of chunk plot-nls" width="400px" />




#### Linear model of transformation
We use were the reciprocal transformation, so C/N versus C.
First we create a the transformed y-variable:

```r
ZINC$Y <- ZINC$C/ZINC$N
```


Fitting a linear model to this data is done with lm():

```r
mod_lm <- lm(Y ~ C, data = ZINC)
plot(ZINC$C, ZINC$Y, ylab = "C/N", xlab = "C")
abline(mod_lm)
```

<img src="figure/plot-lm.png" title="plot of chunk plot-lm" alt="plot of chunk plot-lm" width="400px" />

```r
summary(mod_lm)
```

```
## 
## Call:
## lm(formula = Y ~ C, data = ZINC)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.006926 -0.001708  0.000268  0.003081  0.003706 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.04351    0.00225    19.3  5.3e-08 ***
## C            0.11400    0.00754    15.1  3.6e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0037 on 8 degrees of freedom
## Multiple R-squared: 0.966,	Adjusted R-squared: 0.962 
## F-statistic:  229 on 1 and 8 DF,  p-value: 3.62e-07
```

We get from this K and M as:

* $K = \frac{slope}{intercept} = \frac{0.114}{0.043} = 2.62$
* $M = \frac{1}{slope} = \frac{1}{0.114} = 8.77$

The R^2 is 0.966.


#### Linear model of transformation with weights
Newman used N^4 / C^2 weighting. So first we need to calculate the weights:

```r
ZINC$WGT = ZINC$N^4/ZINC$C^2
```


And fit the linear model with weighting:

```r
mod_wgt <- lm(Y ~ C, data = ZINC, weights = ZINC$WGT)
summary(mod_wgt)
```

```
## 
## Call:
## lm(formula = Y ~ C, data = ZINC, weights = ZINC$WGT)
## 
## Weighted Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.1911 -0.0834  0.0291  0.0580  0.0858 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.04708    0.00199    23.6  1.1e-08 ***
## C            0.10373    0.00568    18.3  8.3e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.105 on 8 degrees of freedom
## Multiple R-squared: 0.977,	Adjusted R-squared: 0.974 
## F-statistic:  333 on 1 and 8 DF,  p-value: 8.32e-08
```

The R^2 is slightly higher: 0.977.

The result for K is:

```r
coef(mod_wgt)[2]/coef(mod_wgt)[1]
```

```
##      C 
## 2.2033
```


and for M:

```r
1/coef(mod_wgt)[2]
```

```
##      C 
## 9.6403
```


#### Are the models appropiate?

We can inspect the residuals of both models:



```r
par(mfrow = c(1, 2))
# lm
plot(mod_lm, which = 1, main = "linear model without weights")
# nls
plot(fitted(mod_nls), residuals(mod_nls), xlab = "fitted", ylab = "Residuals", 
    main = "nonlinear regression")
abline(h = 0, lty = "dotted")
```

<img src="figure/plot-resid.png" title="plot of chunk plot-resid" alt="plot of chunk plot-resid" width="500px" />


The linear model clearly shows an arc-pattern in the residuals - so the data may not follow a linear relationship.
The nonlinear model performs better.



Once again we reproduced the same results as in the book using R :)
Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p101'.


### Quantitative Ecotoxicology, Page 108, Example 3.7, Accumulation

This is example 3.7 on page 108 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about accumulation in mosquitofish (*Gambusia holbrooki*).

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p108.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p108.csv",
ssl.verifypeer = FALSE)
MERCURY <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(MERCURY)
```

```
##   DAY  HG
## 1   0   0
## 2   1 380
## 3   2 540
## 4   3 570
## 5   4 670
## 6   6 780
```


This is pretty much like the previous examples: 

We fit a nonlinear model to our data
.
The model is given in equation 3.42 of the book:

$$C_t = \frac{k_u}{k_e} C_1 (1-e^{-k_e t})$$


```r
plot(MERCURY)
```

<img src="figure/plot_raw.png" title="plot of chunk plot_raw" alt="plot of chunk plot_raw" width="400px" />


We can specify the model as follows:

```r
mod <- nls(HG ~ KU/KE * 0.24 * (1 - exp(-KE * DAY)), data = MERCURY, start = list(KU = 1000, 
    KE = 0.5))
```


This equals to equation 3.42:

* $HG = C_t$
* $KU = k_u$
* $KE = k_e$
* $0.24 = C_1$
* $DAY = t$


Unlike in the book I did not specify bounds here (see the previous posts how to do this).

This results in:

```r
summary(mod)
```

```
## 
## Formula: HG ~ KU/KE * 0.24 * (1 - exp(-KE * DAY))
## 
## Parameters:
##    Estimate Std. Error t value Pr(>|t|)   
## KU 1866.700    241.784    7.72   0.0015 **
## KE    0.589      0.106    5.55   0.0051 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 43.7 on 4 degrees of freedom
## 
## Number of iterations to convergence: 7 
## Achieved convergence tolerance: 2.03e-06
```

So the parameter estimates are:

* $k_e = 0.589 \pm 0.106$
* $k_u = 1866.7 \pm 241.784$

The BCF is given as $BCF = \frac{k_u}{k_e} = 3171.4$

```r
BCF = coef(mod)[1]/coef(mod)[2]
BCF
```

```
##     KU 
## 3171.4
```


From this we can predict the fish concentration as $$C_{fish}=BCF \cdot C_1=761.14$$

```r
BCF * 0.24
```

```
##     KU 
## 761.14
```


Finally we plot the data and our model:

```r
DAY_pred <- seq(0, 6, by = 0.1)
# Raw data
plot(MERCURY)
# add model
lines(DAY_pred, predict(mod, newdata = data.frame(DAY = DAY_pred)))
# add model-equation
text(3, 100, bquote(HG == .(BCF * 0.24) %.% (1 - exp(-.(coef(mod)[2]) %.% DAY))))
```

<img src="figure/plot_model.png" title="plot of chunk plot_model" alt="plot of chunk plot_model" width="400px" />



Once again we reproduced the results as in the book using R :)
The differences for BCF and $C_{fish}$ are due to rounding errors.


Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p108'.
### Quantitative Ecotoxicology, Page 109, Example 3.8, Bioaccumulation

This is example 3.8 on page 109 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about accumulation and elimination of bromophos from water in a guppy (*Poecilia reticulata*).

There are two data files for this example - one for the [accumulation](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_accum.csv) and on for the [elimination](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_elimin.csv).


### Accumulation
First we will look at the accumulation phase:

```r
require(RCurl)
# Accumulation
url_accum <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_accum.csv",
ssl.verifypeer = FALSE)
ACCUM <- read.table(text = url_accum, header = TRUE, sep = ";")
```


```r
head(ACCUM)
```

```
##   HOUR BRPHOS
## 1  0.5   1900
## 2  1.0   3000
## 3  2.0   5200
## 4  4.0   6900
## 5  8.0  24000
## 6 24.0  50000
```


Again we have two columns: One for the time and one for the concentration.


We fit can same model as in [example 3.7](http://edild.github.com/blog/2013/02/24/quant-ecotox-11/) to this data. The uptake $(k_u)$ and elimination $(k_e)$ constants are estimated simultaneously (at the same time):



```r
mod_accum <- nls(BRPHOS ~ KU/KE * 10.5 * (1 - exp(-KE * HOUR)), data = ACCUM, 
    start = list(KU = 100, KE = 0.01))
```

Note that I used different starting values than in the SAS-Code (must be a typo in the book). Also I didn't specify any bounds.

```r
summary(mod_accum)
```

```
## 
## Formula: BRPHOS ~ KU/KE * 10.5 * (1 - exp(-KE * HOUR))
## 
## Parameters:
##     Estimate Std. Error t value Pr(>|t|)    
## KU 344.79786   31.85529   10.82  4.7e-06 ***
## KE   0.00525    0.00103    5.09  0.00094 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 18900 on 8 degrees of freedom
## 
## Number of iterations to convergence: 7 
## Achieved convergence tolerance: 3.78e-06
```



```r
HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1)
# Raw data
plot(ACCUM, main = "Accumulation")
# add model
lines(HOUR_pred, predict(mod_accum, newdata = data.frame(HOUR = HOUR_pred)))
```

<img src="figure/plot_accum_model.png" title="plot of chunk plot_accum_model" alt="plot of chunk plot_accum_model" width="400px" />


So from the accumulation data we estimated the uptake and elimination constants as:

* $k_e = 0.0053 \pm 0.0010$
* $k_u = 344.798 \pm 31.855$




### Sequential estimation
However we could also estimate the elimination constant $(k_e)$ from the elimination phase and then use this estimate for our accumulation data. 

* First estimate $k_e$ from a linear model (linear transformation)
* Plug this estimated $k_e$ into a nonlinear model to estimate $k_u$



```r
# Elimination data
url_elimin <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_elimin.csv")
ELIMIN <- read.table(text = url_elimin, header = TRUE, sep = ";")
```


```r
head(ELIMIN)
```

```
##   HOUR BRPHOS
## 1    0 500000
## 2   12 450000
## 3   24 370000
## 4   48 290000
## 5   72 190000
## 6   96 150000
```

```r
plot(ELIMIN)
```

<img src="figure/plot_elimin_raw.png" title="plot of chunk plot_elimin_raw" alt="plot of chunk plot_elimin_raw" width="400px" />



We will estimate $k_e$ from a linear model like in [previous examples](http://edild.github.com/blog/2013/02/24/quant-ecotox-10/). We could also use nls for this.

First we need to transform the bromophos-concentration to linearize the relationship.

```r
ELIMIN$LBROMO <- log(ELIMIN$BRPHOS)
```


The we can use lm() to fit the linear model:

```r
mod_elimin_lm <- lm(LBROMO ~ HOUR, data = ELIMIN)
summary(mod_elimin_lm)
```

```
## 
## Call:
## lm(formula = LBROMO ~ HOUR, data = ELIMIN)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.09025 -0.03880 -0.00931  0.05900  0.11601 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 13.21262    0.03604   366.6  3.0e-16 ***
## HOUR        -0.01469    0.00025   -58.7  1.1e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0752 on 7 degrees of freedom
## Multiple R-squared: 0.998,	Adjusted R-squared: 0.998 
## F-statistic: 3.44e+03 on 1 and 7 DF,  p-value: 1.09e-10
```


So we get an estimate of $k_e$ as $0.0147 \pm 0.0003$.

This is quite different to the $k_e$ estimated simultaneous from the accumulation data!
Our linear model fits very good (R^2 = 0.998, no pattern in the residuals), so something is strange here...

```r
par(mfrow = c(1, 2))
# plot linearized model
plot(LBROMO ~ HOUR, data = ELIMIN, main = "Data + Model")
# add regression line
abline(mod_elimin_lm)
# plot residuals
plot(fitted(mod_elimin_lm), residuals(mod_elimin_lm), main = "Residuals")
abline(h = 0, lty = "dotted")
```

<img src="figure/elimin_diag.png" title="plot of chunk elimin_diag" alt="plot of chunk elimin_diag" width="400px" />



### Plug $k_e$ from the elimination phase into the accumulation model

Lets take $k_e$ from the elimination phase and plug it into our accumulation model and investigate the differences:


```r
mod_accum2 <- nls(BRPHOS ~ KU/-coef(mod_elimin_lm)[2] * 10.5 * (1 - exp(coef(mod_elimin_lm)[2] * 
    HOUR)), data = ACCUM, start = list(KU = 100))
summary(mod_accum2)
```

```
## 
## Formula: BRPHOS ~ KU/-coef(mod_elimin_lm)[2] * 10.5 * (1 - exp(coef(mod_elimin_lm)[2] * 
##     HOUR))
## 
## Parameters:
##    Estimate Std. Error t value Pr(>|t|)    
## KU    643.9       40.4    15.9  6.7e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 51700 on 9 degrees of freedom
## 
## Number of iterations to convergence: 1 
## Achieved convergence tolerance: 5.27e-09
```


This estimates $k_u = 643.9 \pm 40.4$ which differs greatly from our initial results!
Lets plot this model and the residuals:

```r
par(mfrow = c(1, 2))
HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1)
# Raw data
plot(ACCUM, main = "Accumulation")
# add model
lines(HOUR_pred, predict(mod_accum2, newdata = data.frame(HOUR = HOUR_pred)))
plot(fitted(mod_accum2), residuals(mod_accum2))
```

<img src="figure/plot_accum_model2.png" title="plot of chunk plot_accum_model2" alt="plot of chunk plot_accum_model2" width="400px" />



The residuals show a clear curve pattern. But we could also look at the residual sum of squares and the AIC to see which model fit better to the accumulation data:


```r
# Residual sum of squares
mod_accum$m$deviance()
```

```
## [1] 2870355583
```

```r
mod_accum2$m$deviance()
```

```
## [1] 24088396565
```

```r
# AIC
AIC(mod_accum)
```

```
## [1] 229.13
```

```r
AIC(mod_accum2)
```

```
## [1] 248.4
```


So the first model seem to better fit to the data. However see the discussion in the book for this example!

Once again we reproduced the results as in the book using R :)

Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p109'.

Quantitative Ecotoxicology, Page 147, Example 4.3, LC50

This is about example 4.3 on page 147 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647). This example is about Calculations of $LC_{50}$ values.
In this post I won't reproduce the SAS-Code since I do not have any experience with SAS PROC PROBIT and I do not fully understand whats happening there.

Instead I will fit Dose-Response-Models using the `drc`-package to this data.


Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p147.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p147.csv",
ssl.verifypeer = FALSE)
salt <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(salt)
```

```
##   DEAD TOTAL CONC
## 1   16    76 10.3
## 2   22    79 10.8
## 3   40    77 11.6
## 4   69    76 13.2
## 5   78    78 15.8
## 6   77    77 20.1
```


So the data consists of number of dead animals (DEAD) from all animals (TOTAL) exposed to a concentration (CONC).
First we create a new column with the proportion of dead animals:


```r
salt$prop <- salt$DEAD/salt$TOTAL
```


Lets have a look at the raw data (note that I use a logarithmic scale for the x-axis):

```r
plot(salt$CONC, salt$prop, xlab = "Concentration", ylab = "Proportion dead", 
    log = "x")
```

<img src="figure/p147_plot_raw.png" title="plot of chunk p147_plot_raw" alt="plot of chunk p147_plot_raw" width="400px" />



I will use the drc-package of Christian Ritz and Jens Strebig to fit dose-response-curves to this data. The main function of this package is `drm`:


Here I fit a two-parameter log-logistic model to the data (see <a href="http://dx.doi.org/10.1002/etc.7">Ritz (2010)</a> for a review of dose-response-models used in ecotoxicology):

```r
require(drc)
mod <- drm(prop ~ CONC, data = salt, fct = LL.2())
```


So the usage is similar to `lm()` or `nls()`, except the `fct` argument. This argument defines the model that is fitted to the data.

We can compare this model with other models using the AIC (the smaller the better). 

Here I compare the 2-parameter log-logistic model with a two-parameter Weibull and a 2-parameter Gompertz model. 

`drc` has for this purpose the `mselect()` function that shows some diagnostics for every model:

* likelihood (the higher the better; however use this only when your models are nested)
* AIC (the lower the better)
* residual variance (the lower the better)


```r
mselect(mod, fctList = list(W1.2(), G.2()))
```

```
##      logLik      IC Lack of fit    Res var
## LL.2 14.613 -23.226          NA 0.00067322
## G.2  11.903 -17.806          NA 0.00166155
## W1.2 10.718 -15.437          NA 0.00246584
```


The LL.2-model has the lowest AIC so I will keep this. 
Lets see how the model looks like:

```r
# raw data
plot(prop ~ CONC, data = salt, xlim = c(9, 21), ylim = c(0, 1), ylab = "Proportion dead", 
    xlab = "Concentration")

conc_pred <- seq(9, 21, 0.1)
lines(conc_pred, predict(mod, newdata = data.frame(CONC = conc_pred)))
```

<img src="figure/p147_plot_mod.png" title="plot of chunk p147_plot_mod" alt="plot of chunk p147_plot_mod" width="400px" />


We can get the $LC_{50}$ with confidence interval from the model using the `ED()` function:

```r
ED(mod, 50, interval = "delta")
```

```
## 
## Estimated effective doses
## (Delta method-based confidence interval(s))
## 
##      Estimate Std. Error   Lower Upper
## 1:50  11.4768     0.0575 11.3171  11.6
```



#### References

<p>Ritz C (2010).
&ldquo;Toward A Unified Approach to Dose-Response Modeling in Ecotoxicology.&rdquo;
<EM>Environmental Toxicology And Chemistry</EM>, <B>29</B>.
ISSN 07307268, <a href="http://dx.doi.org/10.1002/etc.7">http://dx.doi.org/10.1002/etc.7</a>.


-----------------------------------
Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p147'.
### Quantitative Ecotoxicology, page 33, example 2.1, Winsorization:

Get the data (Sulfate Concentrations from Savannah River (South Carolina) in mg / L)) from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p33.csv) and read it into R:





```r
ALL <- read.table("p33.csv", header = TRUE, sep = ";")
```


So we have a data.frame with one variable and 21 observations:

```r
str(ALL)
```

```
## 'data.frame':	21 obs. of  1 variable:
##  $ SO4: num  1.3 2.3 2.6 3.3 3.5 3.5 3.6 4 4.1 4.5 ...
```

```r
ALL$SO4
```

```
##  [1] 1.3 2.3 2.6 3.3 3.5 3.5 3.6 4.0 4.1 4.5 5.2 5.6 5.7 6.1 6.2 6.5 6.9
## [18] 7.1 7.7 7.9 9.9
```



Winsorization replaces extreme data values with less extreme values. I have written a small function to run the winsorisation:

```r
winsori <- function(x, width = 2) {
    # check if sorted
    if (is.unsorted(x)) 
        stop("Values must be sorted!")
    # get number of observations
    n <- length(x)
    # Replace lowest
    x[1:width] <- x[width + 1]
    # replace highest
    x[(n - width + 1):n] <- x[(n - width)]
    x
}
```


The function takes a ordered vector and replaces the 2 highest and 2 lowest values (can be changed by the 'width'-Argument by their neighbors.

We can apply this function to our data and safe it as new column:

```r
ALL$SO4_win <- winsori(ALL$SO4)
# display the first and 5 last rows
ALL[c(1:5, 17:21), ]
```

```
##    SO4 SO4_win
## 1  1.3     2.6
## 2  2.3     2.6
## 3  2.6     2.6
## 4  3.3     3.3
## 5  3.5     3.5
## 17 6.9     6.9
## 18 7.1     7.1
## 19 7.7     7.7
## 20 7.9     7.7
## 21 9.9     7.7
```


Worked as expected.
The Winsorized mean and standard-deviation is:

```r
# mean
mean(ALL$SO4_win)
```

```
## [1] 5.081
```

```r
# standard deviation
sd(ALL$SO4_win)
```

```
## [1] 1.792
```


For the Winsorized Standard Deviation we need again a homemade function:

```r
sw <- function(x, width = 2) {
    n <- length(x)
    sd(x) * (n - 1)/(n - 2 * width - 1)
}
sw(ALL$SO4_win)
```

```
## [1] 2.24
```


And lastly we calculate the mean for the trimmed data (remove two observation from each tail):

```r
mean(ALL$SO4, trim = 2/21)
```

```
## [1] 5.065
```


Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p33'.
### Quantitative Ecotoxicology, page 35, Robust Regression on Order Statistics:

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/S35.csv) and read it into R:




```r
SO4 <- read.table("S35.csv", header = TRUE, sep = ";")
```



First we need to convert the vector indicating if an observation is censored to TRUE/FALSE:
I store it in a new colum called 'rem2' (you could also overwrite df$rem):

```r
SO4$rem2 <- ifelse(SO4$rem == "<", TRUE, FALSE)
SO4
```

```
##    value rem  rem2
## 1    2.5   <  TRUE
## 2    2.5   <  TRUE
## 3    2.6   X FALSE
## 4    3.3   X FALSE
## 5    3.5   X FALSE
## 6    3.5   X FALSE
## 7    3.6   X FALSE
## 8    4.0   X FALSE
## 9    4.1   X FALSE
## 10   4.5   X FALSE
## 11   5.2   X FALSE
## 12   5.6   X FALSE
## 13   5.7   X FALSE
## 14   6.1   X FALSE
## 15   6.2   X FALSE
## 16   6.5   X FALSE
## 17   6.9   X FALSE
## 18   7.1   X FALSE
## 19   7.7   X FALSE
## 20   7.9   X FALSE
## 21   9.9   X FALSE
```


Then we can run the Robust Regression on Order Statistics with the ros() function from the NADA package:

```r
require(NADA)
rs <- ros(SO4$value, SO4$rem2)
print(rs)
```

```
##      n  n.cen median   mean     sd 
## 21.000  2.000  5.200  5.158  2.071
```


Which gives the same mean and standard deviation as the SAS-Makro (5.16 and 2.07).

Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p35'.

### Quantitative Ecotoxicology, page 39, example 2.3, Kaplan–Meier estimates

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/S39.csv) and read it into R:




```r
SULFATE <- read.table("S39.csv", header = TRUE, sep = ";")
```


Convert left to right censored data:

```r
SULFATE$FLIP <- abs(SULFATE$SO4 - 8)
SULFATE
```

```
##    SO4 FLAG FLIP
## 1  7.9    1  0.1
## 2  7.7    1  0.3
## 3  7.1    1  0.9
## 4  6.9    1  1.1
## 5  6.5    1  1.5
## 6  6.2    1  1.8
## 7  6.1    1  1.9
## 8  5.7    1  2.3
## 9  5.6    1  2.4
## 10 5.2    1  2.8
## 11 4.5    1  3.5
## 12 4.1    1  3.9
## 13 4.0    1  4.0
## 14 3.6    1  4.4
## 15 3.5    1  4.5
## 16 3.5    1  4.5
## 17 3.3    1  4.7
## 18 2.6    1  5.4
## 19 2.5    0  5.5
## 20 2.5    0  5.5
```


The Kaplan-Meier estimates can be calculated using survfit() from the survival package:

```r
require(survival)
fit <- survfit(Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type = "plain")
fit
```

```
## Call: survfit(formula = Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type = "plain")
## 
## records   n.max n.start  events  median 0.95LCL 0.95UCL 
##   20.00   20.00   20.00   18.00    3.15    1.80    4.50
```


I set conf.type="plain" to be concordant with 'CONFTYPE=LINEAR' from SAS.

The median of 3.15, 95% CI [1.8, 4.5] is the same as with SAS.

Finally a quick plot:

```r
plot(fit)
```

![plot of chunk p39](figure/p39.png) 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p39'.


### Quantitative Ecotoxicology, page 42, example 2.4, Wilcoxon rank sum test

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p42.csv) and read it into R:




```r
SULFATE <- read.table("p42.csv", header = TRUE, sep = ";")
```


Lets first have a look at the data via a violin plot:

```r
require(ggplot2)
ggplot(SULFATE, aes(x = SITE, y = SO4)) + geom_violin()
```

![plot of chunk p43](figure/p43.png) 



It is quite easy to perform a wilcoxon-test with the function wilcox.test:

```r
wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE)
```

```
## Warning: cannot compute exact p-value with ties
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  SO4 by SITE 
## W = 330.5, p-value = 0.00563
## alternative hypothesis: true location shift is not equal to 0
```

It works with the usual formula-notation, additional I specified the continuity correction.
For a one-sided test we can specify the argument 'alternative':

```r
wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE, alternative = "greater")
```

```
## Warning: cannot compute exact p-value with ties
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  SO4 by SITE 
## W = 330.5, p-value = 0.002815
## alternative hypothesis: true location shift is greater than 0
```


The p-values are the same as with SAS, however we get a warning since we have ties in our data (I don't know how SAS handles this):
> cannot compute exact p-value with ties

If we want to compute exact p-values in the presence of ties, we could use wilcox_test() from the coin package: 


```r
require(coin)
wilcox_test(SO4 ~ SITE, SULFATE, distribution = "exact", conf.int = TRUE)
```

```
## 
## 	Exact Wilcoxon Mann-Whitney Rank Sum Test
## 
## data:  SO4 by SITE (A, B) 
## Z = 2.781, p-value = 0.004727
## alternative hypothesis: true mu is not equal to 0 
## 95 percent confidence interval:
##  0.4 2.8 
## sample estimates:
## difference in location 
##                    1.5
```

Here I also specified to output the 95%-Confidence-Interval via the 'conf.int argument.

Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p42'.
### Quantitative Ecotoxicology, page 45, example 2.5, Gehan-Test

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p45.csv) and read it into R:




```r
CADMIUM <- read.table("p45.csv", header = TRUE, sep = ";")
```


'Flip' the data:

```r
CADMIUM$FLIP <- abs(CADMIUM$CD - 100)
CADMIUM
```

```
##      CD SITE FLAG FLIP
## 1  81.3    A    1 18.7
## 2   4.9    A    1 95.1
## 3   4.6    A    1 95.4
## 4   3.5    A    1 96.5
## 5   3.4    A    1 96.6
## 6   3.0    A    1 97.0
## 7   2.9    A    1 97.1
## 8   1.4    B    1 98.6
## 9   0.8    B    1 99.2
## 10  0.7    B    1 99.3
## 11  0.6    A    1 99.4
## 12  0.6    A    1 99.4
## 13  0.6    B    0 99.4
## 14  0.4    B    1 99.6
## 15  0.4    B    1 99.6
## 16  0.4    B    1 99.6
## 17  0.4    B    0 99.6
## 18  0.3    B    0 99.7
## 19  0.2    A    0 99.8
```


And test for differences using survdiff from the survival package:

```r
require(survival)
# log-rank
fit <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
fit
```

```
## Call:
## survdiff(formula = Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
## 
##         N Observed Expected (O-E)^2/E (O-E)^2/V
## SITE=A 10        9     4.99      3.23      5.53
## SITE=B  9        6    10.01      1.61      5.53
## 
##  Chisq= 5.5  on 1 degrees of freedom, p= 0.0187
```

```r
# Peto & Peto modification of the Gehan-Wilcoxon test
fit2 <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
fit2
```

```
## Call:
## survdiff(formula = Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
## 
##         N Observed Expected (O-E)^2/E (O-E)^2/V
## SITE=A 10     6.84     3.55      3.05      7.02
## SITE=B  9     2.84     6.13      1.76      7.02
## 
##  Chisq= 7  on 1 degrees of freedom, p= 0.00808
```




Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p45'.
### Quantitative Ecotoxicology, page 85, example 3.2, Nonlinear Regression

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p85.csv) and read it into R:




```r
OYSTERZN <- read.table("p85.csv", header = TRUE, sep = ";")
```


```r
head(OYSTERZN)
```

```
##   DAY ZINC
## 1   1  700
## 2   1  695
## 3   1  675
## 4   1  630
## 5   1  606
## 6   1  540
```


```r
plot(ZINC ~ DAY, data = OYSTERZN)
```

![plot of chunk p85_raw](figure/p85_raw.png) 


First we fit a **nonlinear Regression without weighting**.

```r
mod_noweight <- nls(ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY), data = OYSTERZN, 
    start = list(KE = 0.004, INITACT = 500))
```


So we fit the model

![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_t%20=%20C_0%20e^{-%28k_{e1}%2Bk_{e2}%29%20t})

to our data.

In the R formula **ZINC** corresponds to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_t) , 
**INITACT** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_0), 
**KE** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e1}), 
**0.00283** is the decay rate constant for 65-Zn ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e2})
and **DAY** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=t).


Were are going to estimate **KE** and **INITACT** and also supplied some start-values for the algorithm.

We can look a the summary to get the estimates and standard error:

```r
summary(mod_noweight)
```

```
## 
## Formula: ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY)
## 
## Parameters:
##         Estimate Std. Error t value Pr(>|t|)    
## KE      2.68e-03   6.68e-04    4.01  0.00015 ***
## INITACT 4.65e+02   2.04e+01   22.86  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 98.8 on 71 degrees of freedom
## 
## Number of iterations to convergence: 3 
## Achieved convergence tolerance: 5.57e-06
```





The resulting estimates of ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e1}) and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_0) are `0.0027` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `6.6841 &times; 10<sup>-4</sup>` and  `465` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `20`.

We can investigate the residuals, which show a clear pattern:


```r
res_mod_noweight <- resid(mod_noweight)
plot(OYSTERZN$DAY, res_mod_noweight)
abline(h = 0)
```

![plot of chunk p85_residuals_nls](figure/p85_residuals_nls.png) 


Secondly, we run a **nonlinear regression with day-squared weighting**:

We use day^2 as weights and add there a column to our data:

```r
OYSTERZN$WNLIN <- OYSTERZN$DAY^2
```


We run again nls, but now we supply this new column as weights:

```r
mod_weight <- nls(ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY), data = OYSTERZN, 
    weights = OYSTERZN$WNLIN, start = list(KE = 0.004, INITACT = 500))
summary(mod_weight)
```

```
## 
## Formula: ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY)
## 
## Parameters:
##         Estimate Std. Error t value Pr(>|t|)    
## KE      2.44e-03   3.47e-04    7.03  1.1e-09 ***
## INITACT 4.55e+02   3.82e+01   11.89  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 7550 on 71 degrees of freedom
## 
## Number of iterations to convergence: 4 
## Achieved convergence tolerance: 5.42e-07
```






The estimates (`0.0024` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `3.4707 &times; 10<sup>-4</sup>` and  `455` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `38` are quite similar to the non weighted regression.

We could plot the two models and the data:


```r
# extract the fitted values
fit_mod_noweight <- fitted(mod_noweight)
fit_mod_weight <- fitted(mod_weight)
# plot data
plot(ZINC ~ DAY, data = OYSTERZN)
# add fitted values
lines(OYSTERZN$DAY, fit_mod_noweight)
lines(OYSTERZN$DAY, fit_mod_weight, lty = "dashed")
# add legend
legend("topright", legend = c("nonweighted", "weighted"), lty = c("solid", "dashed"))
```

![plot of chunk p85_nls_fitted](figure/p85_nls_fitted.png) 



Finally we can also fit a **linear model** to the transformed Zinc-Concentrations:

First we ln-transform the concentrations:

```r
OYSTERZN$LZINC <- log(OYSTERZN$ZINC)
```


We see that the data has now linear trend:

```r
plot(LZINC ~ DAY, OYSTERZN)
```

![plot of chunk p85_raw2](figure/p85_raw2.png) 


And fit a linear regression:

```r
mod_lm <- lm(LZINC ~ DAY, data = OYSTERZN)
sum_lm <- summary(mod_lm)
sum_lm
```

```
## 
## Call:
## lm(formula = LZINC ~ DAY, data = OYSTERZN)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.8974 -0.2448 -0.0709  0.2958  0.5715 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  6.073833   0.056834   106.9   <2e-16 ***
## DAY         -0.005314   0.000243   -21.9   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.354 on 71 degrees of freedom
## Multiple R-squared: 0.871,	Adjusted R-squared: 0.869 
## F-statistic:  478 on 1 and 71 DF,  p-value: <2e-16
```


which is fitting the model
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=ln%28Zn%29%20=%20a%20*%20day%2Bintercept) with a = `-0.0053` and intercept = `6.07`

Now plot data and model, as well as the residuals:

```r
# fitted values
fit_mod_lm <- fitted(mod_lm)

# data + fitted
plot(LZINC ~ DAY, OYSTERZN)
lines(OYSTERZN$DAY, fit_mod_lm)
```

![plot of chunk p85_lm_fitted](figure/p85_lm_fitted1.png) 

```r

# residuals
plot(mod_lm, which = 1)
```

![plot of chunk p85_lm_fitted](figure/p85_lm_fitted2.png) 


The mean square error can be calculated from the summary:

```r
# MSE
sum_lm$sigma^2
```

```
## [1] 0.1252
```

From which we can get an unbiased estimate of $C_0$:

```r
# unbiased estimate for C_0: exp(MSE/2) * exp(Intercept)
exp(sum_lm$sigma^2/2) * exp(sum_lm$coefficients[1, 1])
```

```
## [1] 462.4
```

where 

```r
sum_lm$coefficients[1, 1]
```

extracts the intercept from the summary.

The estimated k in the summary output is `-0.0053` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `2.4298 &times; 10<sup>-4</sup>`, and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_e%20=%20k%20-%20decayrate%20=%200.00531%20-%200.00283%20=%200.00248).

This result is similar to the weighted and non weighted nonlinear regression.
Again we have the same results as with SAS :) [Small deviations may be due to rounding error]




Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p85'.
### Quantitative Ecotoxicology, page 85, example 3.3, Backstripping

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p89.csv) and read it into R:


```r
MERCURY <- read.table("p89.csv", header = TRUE, sep = ";")
```




```r
head(MERCURY)
```

```
##   DAY   PCI
## 1   3 27400
## 2   6 17601
## 3  10 12950
## 4  14 11157
## 5  22  9410
## 6  31  9150
```


```r
plot(PCI ~ DAY, data = MERCURY)
abline(v = 22)
```

![plot of chunk p89_raw](figure/p89_raw.png) 


We can identify two compartments on this plot: A slow one after day 22 and a fast one before day 22.

First we estimate ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_B) and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_B) for the slow compartment, using linear regression of ln-transformed activity against day and predict from this slow compartment the activity over the whole period:


```r
# ln-transformation
MERCURY$LPCI <- log(MERCURY$PCI)
# fit linear model for day 31 to 94
mod_slow <- lm(LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
sum_mod_slow <- summary(mod_slow)
sum_mod_slow
```

```
## 
## Call:
## lm(formula = LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2502 -0.0657  0.0605  0.0822  0.1386 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  9.43096    0.13987   67.42  2.6e-12 ***
## DAY         -0.01240    0.00213   -5.82   0.0004 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.135 on 8 degrees of freedom
## Multiple R-squared: 0.809,	Adjusted R-squared: 0.785 
## F-statistic: 33.9 on 1 and 8 DF,  p-value: 0.000395
```

```r
exp(coef(mod_slow)[1])
```

```
## (Intercept) 
##       12468
```

So this gives us the model ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C = 12468 e^{-0.0124 * Day}) for the slow component.

We do a bias correction and predict the activity for the whole data-range:

```r
# bias-correction
corr_mod_slow <- exp(sum_mod_slow$sigma^2/2)
# add bias corrected predicts to data.frame predict takes the whole
# data.frame es newdata, so we get predicts for every day.
MERCURY$PRED <- exp(predict(mod_slow, newdata = MERCURY)) * corr_mod_slow
# save C_B and k_B as objects (used later...)
CB <- exp(coef(mod_slow)[1]) * corr_mod_slow
KCB <- abs(coef(mod_slow)[2])
```


The residuals from these predictions for day 3 to 22 are associated with the fast compartment.
And we fit a linear regression to the ln-transformed residuals for the fast compartment.


```r
plot(LPCI ~ DAY, data = MERCURY)
points(MERCURY$DAY[1:4], MERCURY$LPCI[1:4], pch = 16)
abline(a = log(CB), b = -KCB)
for (i in 1:4) {
    lines(c(MERCURY$DAY[i], MERCURY$DAY[i]), c(log(MERCURY$PRED[i]), MERCURY$LPCI[i]), 
        lwd = 2)
}
```

![plot of chunk p89_residuals](figure/p89_residuals.png) 



```r
# extract residuals
MERCURY$ERROR <- MERCURY$PCI - MERCURY$PRED
MERCURY$ERROR[1:4]
```

```
## [1] 15276.20  5920.03  1834.41   579.42
```

```r
# fit linear model to ln(residuals) for Day 3 to 22
mod_fast <- lm(log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY < 22, ])
sum_mod_fast <- summary(mod_fast)
sum_mod_fast
```

```
## 
## Call:
## lm(formula = log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY < 22, 
##     ])
## 
## Residuals:
##       1       2       3       4 
##  0.0278 -0.0304 -0.0157  0.0183 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 10.49602    0.03755   279.5 0.000013 ***
## DAY         -0.29659    0.00407   -72.9  0.00019 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0337 on 2 degrees of freedom
## Multiple R-squared:    1,	Adjusted R-squared: 0.999 
## F-statistic: 5.32e+03 on 1 and 2 DF,  p-value: 0.000188
```

```r
exp(coef(mod_fast)[1])
```

```
## (Intercept) 
##       36171
```


So the model for the fast component is: ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C = 36171 e^{-0.297 * Day})


```r
# bias correction
corr_mod_fast <- exp(sum_mod_fast$sigma^2/2)
# save C_A and k_A as objects
CA <- exp(coef(mod_fast)[1]) * corr_mod_fast
KCA <- abs(coef(mod_fast)[2])
```



Now we have two models: one for the fast component and one for the slow component, and we can make a plot similar to Figure 8.1 in Newman and Clements (2008, pp. 119–120). 


```r
plot(LPCI ~ DAY, data = MERCURY)
abline(mod_slow)
abline(mod_fast, lty = "dotted")
legend("topright", c("slow", "backstripped-fast"), lty = c("solid", "dashed"), 
    cex = 0.8)
```

![plot of chunk p89_backstripped](figure/p89_backstripped.png) 

```r
# Estimates
c(CA, KCA, CB, KCB)
```

```
## (Intercept)         DAY (Intercept)         DAY 
##  3.6192e+04  2.9659e-01  1.2583e+04  1.2403e-02
```


We can use this estimates as start-values to fit a non-linear Model to the data (therefore we stored them into objects).

We want to fit the following model:
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C=C_A * e^{-k_A*Day}%2BC_B e^{-k_B*Day})


```r
nls_mod1 <- nls(PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY), 
                data = MERCURY, 
                algorithm = "port",    # to use the bonds
                start = list(KCB = KCB, KCA = KCA, CB = CB, CA = CA),
                lower = c(0, 0, 5000, 20000), 
                upper = c(1, 1, 20000, 45000))
sum_nls_mod1 <- summary(nls_mod1)
sum_nls_mod1
```

```
## 
## Formula: PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY)
## 
## Parameters:
##     Estimate Std. Error t value Pr(>|t|)    
## KCB 1.19e-02   1.41e-03    8.45  3.9e-06 ***
## KCA 2.97e-01   4.47e-02    6.66  3.6e-05 ***
## CB  1.22e+04   8.58e+02   14.27  1.9e-08 ***
## CA  3.79e+04   4.82e+03    7.86  7.7e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 701 on 11 degrees of freedom
## 
## Algorithm "port", convergence message: relative convergence (4)
```

* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_A) is estimated as `37909` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `4824`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_B) is estimated as `12245` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `858`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_A) is estimated as `0.297` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `0.045`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_A) is estimated as `0.012` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `0.001`


And finally we plot data and model.

```r
plot(PCI ~ DAY, data = MERCURY, type = "n")
points(MERCURY$DAY, MERCURY$PCI, pch = ifelse(MERCURY$DAY <= 22, 16, 17))
# smooth line
pred_nls_mod1 <- predict(nls_mod1, newdata = data.frame(DAY = seq(0, 100, 1)))
lines(seq(0, 100, 1), pred_nls_mod1)
legend("topright", c("Fast", "slow"), pch = c(16, 17))
```

![plot of chunk p89_nls](figure/p89_nls.png) 



Again we get nearly the same results with R, except for some differences in the linear models.

This is probably due to the bias-correction in slow-component-model.
We have a MSE of

```r
sum_mod_slow$sigma^2
```

```
## [1] 0.018349
```

which is identical to the book. From the previous example, the bias can be estimated as
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=e^{MSE/2}):

```r
exp(sum_mod_slow$sigma^2/2)
```

```
## [1] 1.0092
```

which is different to the book (1.002).


I have no SAS at hand, so I cannot check this with SAS. However let me know if there is an error
in my calculations.


**Refs**

> Newman, Michael C., and William Henry Clements. Ecotoxicology: A Comprehensive Treatment. Boca Raton: Taylor /& Francis, 2008. 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p89'.
### Quantitative Ecotoxicology, Page 94, Example 3.5


Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p94.csv) and read it into R:


```r
LEAD <- read.table("p94.csv", header = TRUE, sep = ";")
```




```r
head(LEAD)
```

```
##    DAY LEAD
## 1 0.16 41.0
## 2 0.16 31.0
## 3 0.16 25.3
## 4 1.00 30.5
## 5 1.00 22.7
## 6 1.00 22.0
```



As always we first take a look at the data:

```r
plot(LEAD ~ DAY, LEAD)
```

![plot of chunk p94_raw](figure/p94_raw.png) 


A simple power model may fit the data:

$$C_t = C_1~t^{−P}$$

We could fit such model as in example 3.3 via Nonlinear Least Squares or we could try to linearize the relationship by a ln-transform  of both DAY and LEAD:


```r
LEAD$LLEAD <- log(LEAD$LEAD)
LEAD$LDAY <- log(LEAD$DAY)
plot(LLEAD ~ LDAY, LEAD)
```

![plot of chunk p94_linear](figure/p94_linear.png) 


Now we can us lm() to estimate the coefficients and check our model:


```r
# fit model
mod <- lm(LLEAD ~ LDAY, data = LEAD)
```


The residuals show no pattern:

```r
plot(mod, which = 1)
```

![plot of chunk p94_residuals](figure/p94_residuals.png) 


From the model-output:

```r
mod_sum <- summary(mod)
mod_sum
```

```
## 
## Call:
## lm(formula = LLEAD ~ LDAY, data = LEAD)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.4568 -0.1789  0.0372  0.1689  0.4169 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   3.0008     0.0641   46.80  < 2e-16 ***
## LDAY         -0.2715     0.0313   -8.67  1.5e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.238 on 22 degrees of freedom
## Multiple R-squared: 0.773,	Adjusted R-squared: 0.763 
## F-statistic: 75.1 on 1 and 22 DF,  p-value: 1.53e-08
```



We see that out fitted model hast the formula:
$$Ln(LEAD) = 3.0008 - 0.272 ln(DAY)$$
with an R-squared of 0.77 and is statistically significant. The standard errors for the two parameters are 0.064 and 0.031.

So our backtransformed model would be:
$$ LEAD = exp(3.0008)~Day^{-0.272} = 20.68~Day^{-0.272}$$

Finally we can also plot our model:

```r
plot(LLEAD ~ LDAY, LEAD)
abline(mod)
```

![plot of chunk p94_model](figure/p94_model.png) 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p94'.

#!/bin/bash 

cat * >> all.md
### Quantitative Ecotoxicology with R

There is a fresh book about [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647):
> Newman, Michael C. (2012). Quantitative Ecotoxicology, Second Edition. CRC Press


Unfortunately all examples are written in SAS, to which I have no access (5.530€ for a licence of SAS Analytics Pro).

Since I am a poor student who likes R and ecotoxicology, I will try to reproduce the examples given in the book.


All code and data are available at my [github-repository](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology).
### Quantitative Ecotoxicology, Page 101, Example 3.6, Langmuir

This is example 3.6 on page 101 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about adsorption and how to fit an adsorption model to data.

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p101.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p101.csv",
ssl.verifypeer = FALSE)
ZINC <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(ZINC)
```

```
##      N     C
## 1 0.75 0.030
## 2 1.40 0.069
## 3 1.95 0.118
## 4 2.51 0.166
## 5 3.03 0.217
## 6 3.53 0.270
```


So we have a data.frame with two columns,
where N = amount adsorbed (mmol) per unit mass (g) and  C = equilibrium concentration in the aqueous phase (mmol/ml).

We want fit a Langmuir Model (Equation 3.28 in the book) to this data. 

The three methods described are:

* Nonlinear Regression
* linear transformation
* linear transformation with weighting



#### Nonlinear Regression

```r
mod_nls <- nls(N ~ (K * C * M)/(1 + K * C), data = ZINC, start = list(K = 3, 
    M = 9), lower = 0, algorithm = "port")
```

This fits the model 

$$ N = \frac{KCM}{1+KC} $$ 

to the data. 

We supplied some starting values and specified the lower bonds for K and M as 0 (bonds can only be used with the port algorithm).

This gives us the estimates for K and M as:

```r
summary(mod_nls)
```

```
## 
## Formula: N ~ (K * C * M)/(1 + K * C)
## 
## Parameters:
##   Estimate Std. Error t value Pr(>|t|)    
## K    2.097      0.188    11.1  3.8e-06 ***
## M    9.899      0.521    19.0  6.1e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0929 on 8 degrees of freedom
## 
## Algorithm "port", convergence message: relative convergence (4)
```


* $K = 2.097 \pm 0.188$
* $M = 9.899 \pm 0.521$

The t and p-values of this output are not of interest for us (tests if the parameters deviate from 0).

We can plot the raw data and the model easily using the predict-function:

```r
plot(ZINC$C, ZINC$N, xlab = "C", ylab = "N")
# generate C-values to predict
x_n <- seq(min(ZINC$C), max(ZINC$C), length.out = 200)
# add predicts to plot
lines(x_n, predict(mod_nls, newdata = data.frame(C = x_n)))
```

<img src="figure/plot-nls.png" title="plot of chunk plot-nls" alt="plot of chunk plot-nls" width="400px" />




#### Linear model of transformation
We use were the reciprocal transformation, so C/N versus C.
First we create a the transformed y-variable:

```r
ZINC$Y <- ZINC$C/ZINC$N
```


Fitting a linear model to this data is done with lm():

```r
mod_lm <- lm(Y ~ C, data = ZINC)
plot(ZINC$C, ZINC$Y, ylab = "C/N", xlab = "C")
abline(mod_lm)
```

<img src="figure/plot-lm.png" title="plot of chunk plot-lm" alt="plot of chunk plot-lm" width="400px" />

```r
summary(mod_lm)
```

```
## 
## Call:
## lm(formula = Y ~ C, data = ZINC)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.006926 -0.001708  0.000268  0.003081  0.003706 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.04351    0.00225    19.3  5.3e-08 ***
## C            0.11400    0.00754    15.1  3.6e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0037 on 8 degrees of freedom
## Multiple R-squared: 0.966,	Adjusted R-squared: 0.962 
## F-statistic:  229 on 1 and 8 DF,  p-value: 3.62e-07
```

We get from this K and M as:

* $K = \frac{slope}{intercept} = \frac{0.114}{0.043} = 2.62$
* $M = \frac{1}{slope} = \frac{1}{0.114} = 8.77$

The R^2 is 0.966.


#### Linear model of transformation with weights
Newman used N^4 / C^2 weighting. So first we need to calculate the weights:

```r
ZINC$WGT = ZINC$N^4/ZINC$C^2
```


And fit the linear model with weighting:

```r
mod_wgt <- lm(Y ~ C, data = ZINC, weights = ZINC$WGT)
summary(mod_wgt)
```

```
## 
## Call:
## lm(formula = Y ~ C, data = ZINC, weights = ZINC$WGT)
## 
## Weighted Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.1911 -0.0834  0.0291  0.0580  0.0858 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.04708    0.00199    23.6  1.1e-08 ***
## C            0.10373    0.00568    18.3  8.3e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.105 on 8 degrees of freedom
## Multiple R-squared: 0.977,	Adjusted R-squared: 0.974 
## F-statistic:  333 on 1 and 8 DF,  p-value: 8.32e-08
```

The R^2 is slightly higher: 0.977.

The result for K is:

```r
coef(mod_wgt)[2]/coef(mod_wgt)[1]
```

```
##      C 
## 2.2033
```


and for M:

```r
1/coef(mod_wgt)[2]
```

```
##      C 
## 9.6403
```


#### Are the models appropiate?

We can inspect the residuals of both models:



```r
par(mfrow = c(1, 2))
# lm
plot(mod_lm, which = 1, main = "linear model without weights")
# nls
plot(fitted(mod_nls), residuals(mod_nls), xlab = "fitted", ylab = "Residuals", 
    main = "nonlinear regression")
abline(h = 0, lty = "dotted")
```

<img src="figure/plot-resid.png" title="plot of chunk plot-resid" alt="plot of chunk plot-resid" width="500px" />


The linear model clearly shows an arc-pattern in the residuals - so the data may not follow a linear relationship.
The nonlinear model performs better.



Once again we reproduced the same results as in the book using R :)
Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p101'.


### Quantitative Ecotoxicology, Page 108, Example 3.7, Accumulation

This is example 3.7 on page 108 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about accumulation in mosquitofish (*Gambusia holbrooki*).

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p108.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p108.csv",
ssl.verifypeer = FALSE)
MERCURY <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(MERCURY)
```

```
##   DAY  HG
## 1   0   0
## 2   1 380
## 3   2 540
## 4   3 570
## 5   4 670
## 6   6 780
```


This is pretty much like the previous examples: 

We fit a nonlinear model to our data
.
The model is given in equation 3.42 of the book:

$$C_t = \frac{k_u}{k_e} C_1 (1-e^{-k_e t})$$


```r
plot(MERCURY)
```

<img src="figure/plot_raw.png" title="plot of chunk plot_raw" alt="plot of chunk plot_raw" width="400px" />


We can specify the model as follows:

```r
mod <- nls(HG ~ KU/KE * 0.24 * (1 - exp(-KE * DAY)), data = MERCURY, start = list(KU = 1000, 
    KE = 0.5))
```


This equals to equation 3.42:

* $HG = C_t$
* $KU = k_u$
* $KE = k_e$
* $0.24 = C_1$
* $DAY = t$


Unlike in the book I did not specify bounds here (see the previous posts how to do this).

This results in:

```r
summary(mod)
```

```
## 
## Formula: HG ~ KU/KE * 0.24 * (1 - exp(-KE * DAY))
## 
## Parameters:
##    Estimate Std. Error t value Pr(>|t|)   
## KU 1866.700    241.784    7.72   0.0015 **
## KE    0.589      0.106    5.55   0.0051 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 43.7 on 4 degrees of freedom
## 
## Number of iterations to convergence: 7 
## Achieved convergence tolerance: 2.03e-06
```

So the parameter estimates are:

* $k_e = 0.589 \pm 0.106$
* $k_u = 1866.7 \pm 241.784$

The BCF is given as $BCF = \frac{k_u}{k_e} = 3171.4$

```r
BCF = coef(mod)[1]/coef(mod)[2]
BCF
```

```
##     KU 
## 3171.4
```


From this we can predict the fish concentration as $$C_{fish}=BCF \cdot C_1=761.14$$

```r
BCF * 0.24
```

```
##     KU 
## 761.14
```


Finally we plot the data and our model:

```r
DAY_pred <- seq(0, 6, by = 0.1)
# Raw data
plot(MERCURY)
# add model
lines(DAY_pred, predict(mod, newdata = data.frame(DAY = DAY_pred)))
# add model-equation
text(3, 100, bquote(HG == .(BCF * 0.24) %.% (1 - exp(-.(coef(mod)[2]) %.% DAY))))
```

<img src="figure/plot_model.png" title="plot of chunk plot_model" alt="plot of chunk plot_model" width="400px" />



Once again we reproduced the results as in the book using R :)
The differences for BCF and $C_{fish}$ are due to rounding errors.


Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p108'.
### Quantitative Ecotoxicology, Page 109, Example 3.8, Bioaccumulation

This is example 3.8 on page 109 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647) - reproduced with R. This example is about accumulation and elimination of bromophos from water in a guppy (*Poecilia reticulata*).

There are two data files for this example - one for the [accumulation](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_accum.csv) and on for the [elimination](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_elimin.csv).


### Accumulation
First we will look at the accumulation phase:

```r
require(RCurl)
# Accumulation
url_accum <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_accum.csv",
ssl.verifypeer = FALSE)
ACCUM <- read.table(text = url_accum, header = TRUE, sep = ";")
```


```r
head(ACCUM)
```

```
##   HOUR BRPHOS
## 1  0.5   1900
## 2  1.0   3000
## 3  2.0   5200
## 4  4.0   6900
## 5  8.0  24000
## 6 24.0  50000
```


Again we have two columns: One for the time and one for the concentration.


We fit can same model as in [example 3.7](http://edild.github.com/blog/2013/02/24/quant-ecotox-11/) to this data. The uptake $(k_u)$ and elimination $(k_e)$ constants are estimated simultaneously (at the same time):



```r
mod_accum <- nls(BRPHOS ~ KU/KE * 10.5 * (1 - exp(-KE * HOUR)), data = ACCUM, 
    start = list(KU = 100, KE = 0.01))
```

Note that I used different starting values than in the SAS-Code (must be a typo in the book). Also I didn't specify any bounds.

```r
summary(mod_accum)
```

```
## 
## Formula: BRPHOS ~ KU/KE * 10.5 * (1 - exp(-KE * HOUR))
## 
## Parameters:
##     Estimate Std. Error t value Pr(>|t|)    
## KU 344.79786   31.85529   10.82  4.7e-06 ***
## KE   0.00525    0.00103    5.09  0.00094 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 18900 on 8 degrees of freedom
## 
## Number of iterations to convergence: 7 
## Achieved convergence tolerance: 3.78e-06
```



```r
HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1)
# Raw data
plot(ACCUM, main = "Accumulation")
# add model
lines(HOUR_pred, predict(mod_accum, newdata = data.frame(HOUR = HOUR_pred)))
```

<img src="figure/plot_accum_model.png" title="plot of chunk plot_accum_model" alt="plot of chunk plot_accum_model" width="400px" />


So from the accumulation data we estimated the uptake and elimination constants as:

* $k_e = 0.0053 \pm 0.0010$
* $k_u = 344.798 \pm 31.855$




### Sequential estimation
However we could also estimate the elimination constant $(k_e)$ from the elimination phase and then use this estimate for our accumulation data. 

* First estimate $k_e$ from a linear model (linear transformation)
* Plug this estimated $k_e$ into a nonlinear model to estimate $k_u$



```r
# Elimination data
url_elimin <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p109_elimin.csv")
ELIMIN <- read.table(text = url_elimin, header = TRUE, sep = ";")
```


```r
head(ELIMIN)
```

```
##   HOUR BRPHOS
## 1    0 500000
## 2   12 450000
## 3   24 370000
## 4   48 290000
## 5   72 190000
## 6   96 150000
```

```r
plot(ELIMIN)
```

<img src="figure/plot_elimin_raw.png" title="plot of chunk plot_elimin_raw" alt="plot of chunk plot_elimin_raw" width="400px" />



We will estimate $k_e$ from a linear model like in [previous examples](http://edild.github.com/blog/2013/02/24/quant-ecotox-10/). We could also use nls for this.

First we need to transform the bromophos-concentration to linearize the relationship.

```r
ELIMIN$LBROMO <- log(ELIMIN$BRPHOS)
```


The we can use lm() to fit the linear model:

```r
mod_elimin_lm <- lm(LBROMO ~ HOUR, data = ELIMIN)
summary(mod_elimin_lm)
```

```
## 
## Call:
## lm(formula = LBROMO ~ HOUR, data = ELIMIN)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.09025 -0.03880 -0.00931  0.05900  0.11601 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 13.21262    0.03604   366.6  3.0e-16 ***
## HOUR        -0.01469    0.00025   -58.7  1.1e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0752 on 7 degrees of freedom
## Multiple R-squared: 0.998,	Adjusted R-squared: 0.998 
## F-statistic: 3.44e+03 on 1 and 7 DF,  p-value: 1.09e-10
```


So we get an estimate of $k_e$ as $0.0147 \pm 0.0003$.

This is quite different to the $k_e$ estimated simultaneous from the accumulation data!
Our linear model fits very good (R^2 = 0.998, no pattern in the residuals), so something is strange here...

```r
par(mfrow = c(1, 2))
# plot linearized model
plot(LBROMO ~ HOUR, data = ELIMIN, main = "Data + Model")
# add regression line
abline(mod_elimin_lm)
# plot residuals
plot(fitted(mod_elimin_lm), residuals(mod_elimin_lm), main = "Residuals")
abline(h = 0, lty = "dotted")
```

<img src="figure/elimin_diag.png" title="plot of chunk elimin_diag" alt="plot of chunk elimin_diag" width="400px" />



### Plug $k_e$ from the elimination phase into the accumulation model

Lets take $k_e$ from the elimination phase and plug it into our accumulation model and investigate the differences:


```r
mod_accum2 <- nls(BRPHOS ~ KU/-coef(mod_elimin_lm)[2] * 10.5 * (1 - exp(coef(mod_elimin_lm)[2] * 
    HOUR)), data = ACCUM, start = list(KU = 100))
summary(mod_accum2)
```

```
## 
## Formula: BRPHOS ~ KU/-coef(mod_elimin_lm)[2] * 10.5 * (1 - exp(coef(mod_elimin_lm)[2] * 
##     HOUR))
## 
## Parameters:
##    Estimate Std. Error t value Pr(>|t|)    
## KU    643.9       40.4    15.9  6.7e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 51700 on 9 degrees of freedom
## 
## Number of iterations to convergence: 1 
## Achieved convergence tolerance: 5.27e-09
```


This estimates $k_u = 643.9 \pm 40.4$ which differs greatly from our initial results!
Lets plot this model and the residuals:

```r
par(mfrow = c(1, 2))
HOUR_pred <- seq(min(ACCUM$HOUR), max(ACCUM$HOUR), by = 0.1)
# Raw data
plot(ACCUM, main = "Accumulation")
# add model
lines(HOUR_pred, predict(mod_accum2, newdata = data.frame(HOUR = HOUR_pred)))
plot(fitted(mod_accum2), residuals(mod_accum2))
```

<img src="figure/plot_accum_model2.png" title="plot of chunk plot_accum_model2" alt="plot of chunk plot_accum_model2" width="400px" />



The residuals show a clear curve pattern. But we could also look at the residual sum of squares and the AIC to see which model fit better to the accumulation data:


```r
# Residual sum of squares
mod_accum$m$deviance()
```

```
## [1] 2870355583
```

```r
mod_accum2$m$deviance()
```

```
## [1] 24088396565
```

```r
# AIC
AIC(mod_accum)
```

```
## [1] 229.13
```

```r
AIC(mod_accum2)
```

```
## [1] 248.4
```


So the first model seem to better fit to the data. However see the discussion in the book for this example!

Once again we reproduced the results as in the book using R :)

Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p109'.

Quantitative Ecotoxicology, Page 147, Example 4.3, LC50

This is about example 4.3 on page 147 of [Quantitative Ecotoxicology](http://www.crcpress.com/product/isbn/9781439835647). This example is about Calculations of $LC_{50}$ values.
In this post I won't reproduce the SAS-Code since I do not have any experience with SAS PROC PROBIT and I do not fully understand whats happening there.

Instead I will fit Dose-Response-Models using the `drc`-package to this data.


Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p147.csv) and read it into R:


```r
require(RCurl)
url <- getURL("https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p147.csv",
ssl.verifypeer = FALSE)
salt <- read.table(text = url, header = TRUE, sep = ";")
```


```r
head(salt)
```

```
##   DEAD TOTAL CONC
## 1   16    76 10.3
## 2   22    79 10.8
## 3   40    77 11.6
## 4   69    76 13.2
## 5   78    78 15.8
## 6   77    77 20.1
```


So the data consists of number of dead animals (DEAD) from all animals (TOTAL) exposed to a concentration (CONC).
First we create a new column with the proportion of dead animals:


```r
salt$prop <- salt$DEAD/salt$TOTAL
```


Lets have a look at the raw data (note that I use a logarithmic scale for the x-axis):

```r
plot(salt$CONC, salt$prop, xlab = "Concentration", ylab = "Proportion dead", 
    log = "x")
```

<img src="figure/p147_plot_raw.png" title="plot of chunk p147_plot_raw" alt="plot of chunk p147_plot_raw" width="400px" />



I will use the drc-package of Christian Ritz and Jens Strebig to fit dose-response-curves to this data. The main function of this package is `drm`:


Here I fit a two-parameter log-logistic model to the data (see <a href="http://dx.doi.org/10.1002/etc.7">Ritz (2010)</a> for a review of dose-response-models used in ecotoxicology):

```r
require(drc)
mod <- drm(prop ~ CONC, data = salt, fct = LL.2())
```


So the usage is similar to `lm()` or `nls()`, except the `fct` argument. This argument defines the model that is fitted to the data.

We can compare this model with other models using the AIC (the smaller the better). 

Here I compare the 2-parameter log-logistic model with a two-parameter Weibull and a 2-parameter Gompertz model. 

`drc` has for this purpose the `mselect()` function that shows some diagnostics for every model:

* likelihood (the higher the better; however use this only when your models are nested)
* AIC (the lower the better)
* residual variance (the lower the better)


```r
mselect(mod, fctList = list(W1.2(), G.2()))
```

```
##      logLik      IC Lack of fit    Res var
## LL.2 14.613 -23.226          NA 0.00067322
## G.2  11.903 -17.806          NA 0.00166155
## W1.2 10.718 -15.437          NA 0.00246584
```


The LL.2-model has the lowest AIC so I will keep this. 
Lets see how the model looks like:

```r
# raw data
plot(prop ~ CONC, data = salt, xlim = c(9, 21), ylim = c(0, 1), ylab = "Proportion dead", 
    xlab = "Concentration")

conc_pred <- seq(9, 21, 0.1)
lines(conc_pred, predict(mod, newdata = data.frame(CONC = conc_pred)))
```

<img src="figure/p147_plot_mod.png" title="plot of chunk p147_plot_mod" alt="plot of chunk p147_plot_mod" width="400px" />


We can get the $LC_{50}$ with confidence interval from the model using the `ED()` function:

```r
ED(mod, 50, interval = "delta")
```

```
## 
## Estimated effective doses
## (Delta method-based confidence interval(s))
## 
##      Estimate Std. Error   Lower Upper
## 1:50  11.4768     0.0575 11.3171  11.6
```



#### References

<p>Ritz C (2010).
&ldquo;Toward A Unified Approach to Dose-Response Modeling in Ecotoxicology.&rdquo;
<EM>Environmental Toxicology And Chemistry</EM>, <B>29</B>.
ISSN 07307268, <a href="http://dx.doi.org/10.1002/etc.7">http://dx.doi.org/10.1002/etc.7</a>.


-----------------------------------
Code and data are available on my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p147'.
### Quantitative Ecotoxicology, page 33, example 2.1, Winsorization:

Get the data (Sulfate Concentrations from Savannah River (South Carolina) in mg / L)) from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p33.csv) and read it into R:





```r
ALL <- read.table("p33.csv", header = TRUE, sep = ";")
```


So we have a data.frame with one variable and 21 observations:

```r
str(ALL)
```

```
## 'data.frame':	21 obs. of  1 variable:
##  $ SO4: num  1.3 2.3 2.6 3.3 3.5 3.5 3.6 4 4.1 4.5 ...
```

```r
ALL$SO4
```

```
##  [1] 1.3 2.3 2.6 3.3 3.5 3.5 3.6 4.0 4.1 4.5 5.2 5.6 5.7 6.1 6.2 6.5 6.9
## [18] 7.1 7.7 7.9 9.9
```



Winsorization replaces extreme data values with less extreme values. I have written a small function to run the winsorisation:

```r
winsori <- function(x, width = 2) {
    # check if sorted
    if (is.unsorted(x)) 
        stop("Values must be sorted!")
    # get number of observations
    n <- length(x)
    # Replace lowest
    x[1:width] <- x[width + 1]
    # replace highest
    x[(n - width + 1):n] <- x[(n - width)]
    x
}
```


The function takes a ordered vector and replaces the 2 highest and 2 lowest values (can be changed by the 'width'-Argument by their neighbors.

We can apply this function to our data and safe it as new column:

```r
ALL$SO4_win <- winsori(ALL$SO4)
# display the first and 5 last rows
ALL[c(1:5, 17:21), ]
```

```
##    SO4 SO4_win
## 1  1.3     2.6
## 2  2.3     2.6
## 3  2.6     2.6
## 4  3.3     3.3
## 5  3.5     3.5
## 17 6.9     6.9
## 18 7.1     7.1
## 19 7.7     7.7
## 20 7.9     7.7
## 21 9.9     7.7
```


Worked as expected.
The Winsorized mean and standard-deviation is:

```r
# mean
mean(ALL$SO4_win)
```

```
## [1] 5.081
```

```r
# standard deviation
sd(ALL$SO4_win)
```

```
## [1] 1.792
```


For the Winsorized Standard Deviation we need again a homemade function:

```r
sw <- function(x, width = 2) {
    n <- length(x)
    sd(x) * (n - 1)/(n - 2 * width - 1)
}
sw(ALL$SO4_win)
```

```
## [1] 2.24
```


And lastly we calculate the mean for the trimmed data (remove two observation from each tail):

```r
mean(ALL$SO4, trim = 2/21)
```

```
## [1] 5.065
```


Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p33'.
### Quantitative Ecotoxicology, page 35, Robust Regression on Order Statistics:

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/S35.csv) and read it into R:




```r
SO4 <- read.table("S35.csv", header = TRUE, sep = ";")
```



First we need to convert the vector indicating if an observation is censored to TRUE/FALSE:
I store it in a new colum called 'rem2' (you could also overwrite df$rem):

```r
SO4$rem2 <- ifelse(SO4$rem == "<", TRUE, FALSE)
SO4
```

```
##    value rem  rem2
## 1    2.5   <  TRUE
## 2    2.5   <  TRUE
## 3    2.6   X FALSE
## 4    3.3   X FALSE
## 5    3.5   X FALSE
## 6    3.5   X FALSE
## 7    3.6   X FALSE
## 8    4.0   X FALSE
## 9    4.1   X FALSE
## 10   4.5   X FALSE
## 11   5.2   X FALSE
## 12   5.6   X FALSE
## 13   5.7   X FALSE
## 14   6.1   X FALSE
## 15   6.2   X FALSE
## 16   6.5   X FALSE
## 17   6.9   X FALSE
## 18   7.1   X FALSE
## 19   7.7   X FALSE
## 20   7.9   X FALSE
## 21   9.9   X FALSE
```


Then we can run the Robust Regression on Order Statistics with the ros() function from the NADA package:

```r
require(NADA)
rs <- ros(SO4$value, SO4$rem2)
print(rs)
```

```
##      n  n.cen median   mean     sd 
## 21.000  2.000  5.200  5.158  2.071
```


Which gives the same mean and standard deviation as the SAS-Makro (5.16 and 2.07).

Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p35'.

### Quantitative Ecotoxicology, page 39, example 2.3, Kaplan–Meier estimates

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/S39.csv) and read it into R:




```r
SULFATE <- read.table("S39.csv", header = TRUE, sep = ";")
```


Convert left to right censored data:

```r
SULFATE$FLIP <- abs(SULFATE$SO4 - 8)
SULFATE
```

```
##    SO4 FLAG FLIP
## 1  7.9    1  0.1
## 2  7.7    1  0.3
## 3  7.1    1  0.9
## 4  6.9    1  1.1
## 5  6.5    1  1.5
## 6  6.2    1  1.8
## 7  6.1    1  1.9
## 8  5.7    1  2.3
## 9  5.6    1  2.4
## 10 5.2    1  2.8
## 11 4.5    1  3.5
## 12 4.1    1  3.9
## 13 4.0    1  4.0
## 14 3.6    1  4.4
## 15 3.5    1  4.5
## 16 3.5    1  4.5
## 17 3.3    1  4.7
## 18 2.6    1  5.4
## 19 2.5    0  5.5
## 20 2.5    0  5.5
```


The Kaplan-Meier estimates can be calculated using survfit() from the survival package:

```r
require(survival)
fit <- survfit(Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type = "plain")
fit
```

```
## Call: survfit(formula = Surv(FLIP, FLAG) ~ 1, data = SULFATE, conf.type = "plain")
## 
## records   n.max n.start  events  median 0.95LCL 0.95UCL 
##   20.00   20.00   20.00   18.00    3.15    1.80    4.50
```


I set conf.type="plain" to be concordant with 'CONFTYPE=LINEAR' from SAS.

The median of 3.15, 95% CI [1.8, 4.5] is the same as with SAS.

Finally a quick plot:

```r
plot(fit)
```

![plot of chunk p39](figure/p39.png) 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p39'.


### Quantitative Ecotoxicology, page 42, example 2.4, Wilcoxon rank sum test

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p42.csv) and read it into R:




```r
SULFATE <- read.table("p42.csv", header = TRUE, sep = ";")
```


Lets first have a look at the data via a violin plot:

```r
require(ggplot2)
ggplot(SULFATE, aes(x = SITE, y = SO4)) + geom_violin()
```

![plot of chunk p43](figure/p43.png) 



It is quite easy to perform a wilcoxon-test with the function wilcox.test:

```r
wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE)
```

```
## Warning: cannot compute exact p-value with ties
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  SO4 by SITE 
## W = 330.5, p-value = 0.00563
## alternative hypothesis: true location shift is not equal to 0
```

It works with the usual formula-notation, additional I specified the continuity correction.
For a one-sided test we can specify the argument 'alternative':

```r
wilcox.test(SO4 ~ SITE, data = SULFATE, correct = TRUE, alternative = "greater")
```

```
## Warning: cannot compute exact p-value with ties
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  SO4 by SITE 
## W = 330.5, p-value = 0.002815
## alternative hypothesis: true location shift is greater than 0
```


The p-values are the same as with SAS, however we get a warning since we have ties in our data (I don't know how SAS handles this):
> cannot compute exact p-value with ties

If we want to compute exact p-values in the presence of ties, we could use wilcox_test() from the coin package: 


```r
require(coin)
wilcox_test(SO4 ~ SITE, SULFATE, distribution = "exact", conf.int = TRUE)
```

```
## 
## 	Exact Wilcoxon Mann-Whitney Rank Sum Test
## 
## data:  SO4 by SITE (A, B) 
## Z = 2.781, p-value = 0.004727
## alternative hypothesis: true mu is not equal to 0 
## 95 percent confidence interval:
##  0.4 2.8 
## sample estimates:
## difference in location 
##                    1.5
```

Here I also specified to output the 95%-Confidence-Interval via the 'conf.int argument.

Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p42'.
### Quantitative Ecotoxicology, page 45, example 2.5, Gehan-Test

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p45.csv) and read it into R:




```r
CADMIUM <- read.table("p45.csv", header = TRUE, sep = ";")
```


'Flip' the data:

```r
CADMIUM$FLIP <- abs(CADMIUM$CD - 100)
CADMIUM
```

```
##      CD SITE FLAG FLIP
## 1  81.3    A    1 18.7
## 2   4.9    A    1 95.1
## 3   4.6    A    1 95.4
## 4   3.5    A    1 96.5
## 5   3.4    A    1 96.6
## 6   3.0    A    1 97.0
## 7   2.9    A    1 97.1
## 8   1.4    B    1 98.6
## 9   0.8    B    1 99.2
## 10  0.7    B    1 99.3
## 11  0.6    A    1 99.4
## 12  0.6    A    1 99.4
## 13  0.6    B    0 99.4
## 14  0.4    B    1 99.6
## 15  0.4    B    1 99.6
## 16  0.4    B    1 99.6
## 17  0.4    B    0 99.6
## 18  0.3    B    0 99.7
## 19  0.2    A    0 99.8
```


And test for differences using survdiff from the survival package:

```r
require(survival)
# log-rank
fit <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
fit
```

```
## Call:
## survdiff(formula = Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 0)
## 
##         N Observed Expected (O-E)^2/E (O-E)^2/V
## SITE=A 10        9     4.99      3.23      5.53
## SITE=B  9        6    10.01      1.61      5.53
## 
##  Chisq= 5.5  on 1 degrees of freedom, p= 0.0187
```

```r
# Peto & Peto modification of the Gehan-Wilcoxon test
fit2 <- survdiff(Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
fit2
```

```
## Call:
## survdiff(formula = Surv(FLIP, FLAG) ~ SITE, data = CADMIUM, rho = 1)
## 
##         N Observed Expected (O-E)^2/E (O-E)^2/V
## SITE=A 10     6.84     3.55      3.05      7.02
## SITE=B  9     2.84     6.13      1.76      7.02
## 
##  Chisq= 7  on 1 degrees of freedom, p= 0.00808
```




Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p45'.
### Quantitative Ecotoxicology, page 85, example 3.2, Nonlinear Regression

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p85.csv) and read it into R:




```r
OYSTERZN <- read.table("p85.csv", header = TRUE, sep = ";")
```


```r
head(OYSTERZN)
```

```
##   DAY ZINC
## 1   1  700
## 2   1  695
## 3   1  675
## 4   1  630
## 5   1  606
## 6   1  540
```


```r
plot(ZINC ~ DAY, data = OYSTERZN)
```

![plot of chunk p85_raw](figure/p85_raw.png) 


First we fit a **nonlinear Regression without weighting**.

```r
mod_noweight <- nls(ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY), data = OYSTERZN, 
    start = list(KE = 0.004, INITACT = 500))
```


So we fit the model

![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_t%20=%20C_0%20e^{-%28k_{e1}%2Bk_{e2}%29%20t})

to our data.

In the R formula **ZINC** corresponds to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_t) , 
**INITACT** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_0), 
**KE** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e1}), 
**0.00283** is the decay rate constant for 65-Zn ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e2})
and **DAY** to ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=t).


Were are going to estimate **KE** and **INITACT** and also supplied some start-values for the algorithm.

We can look a the summary to get the estimates and standard error:

```r
summary(mod_noweight)
```

```
## 
## Formula: ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY)
## 
## Parameters:
##         Estimate Std. Error t value Pr(>|t|)    
## KE      2.68e-03   6.68e-04    4.01  0.00015 ***
## INITACT 4.65e+02   2.04e+01   22.86  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 98.8 on 71 degrees of freedom
## 
## Number of iterations to convergence: 3 
## Achieved convergence tolerance: 5.57e-06
```





The resulting estimates of ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_{e1}) and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_0) are `0.0027` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `6.6841 &times; 10<sup>-4</sup>` and  `465` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `20`.

We can investigate the residuals, which show a clear pattern:


```r
res_mod_noweight <- resid(mod_noweight)
plot(OYSTERZN$DAY, res_mod_noweight)
abline(h = 0)
```

![plot of chunk p85_residuals_nls](figure/p85_residuals_nls.png) 


Secondly, we run a **nonlinear regression with day-squared weighting**:

We use day^2 as weights and add there a column to our data:

```r
OYSTERZN$WNLIN <- OYSTERZN$DAY^2
```


We run again nls, but now we supply this new column as weights:

```r
mod_weight <- nls(ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY), data = OYSTERZN, 
    weights = OYSTERZN$WNLIN, start = list(KE = 0.004, INITACT = 500))
summary(mod_weight)
```

```
## 
## Formula: ZINC ~ INITACT * exp((-(KE + 0.00283)) * DAY)
## 
## Parameters:
##         Estimate Std. Error t value Pr(>|t|)    
## KE      2.44e-03   3.47e-04    7.03  1.1e-09 ***
## INITACT 4.55e+02   3.82e+01   11.89  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 7550 on 71 degrees of freedom
## 
## Number of iterations to convergence: 4 
## Achieved convergence tolerance: 5.42e-07
```






The estimates (`0.0024` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `3.4707 &times; 10<sup>-4</sup>` and  `455` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `38` are quite similar to the non weighted regression.

We could plot the two models and the data:


```r
# extract the fitted values
fit_mod_noweight <- fitted(mod_noweight)
fit_mod_weight <- fitted(mod_weight)
# plot data
plot(ZINC ~ DAY, data = OYSTERZN)
# add fitted values
lines(OYSTERZN$DAY, fit_mod_noweight)
lines(OYSTERZN$DAY, fit_mod_weight, lty = "dashed")
# add legend
legend("topright", legend = c("nonweighted", "weighted"), lty = c("solid", "dashed"))
```

![plot of chunk p85_nls_fitted](figure/p85_nls_fitted.png) 



Finally we can also fit a **linear model** to the transformed Zinc-Concentrations:

First we ln-transform the concentrations:

```r
OYSTERZN$LZINC <- log(OYSTERZN$ZINC)
```


We see that the data has now linear trend:

```r
plot(LZINC ~ DAY, OYSTERZN)
```

![plot of chunk p85_raw2](figure/p85_raw2.png) 


And fit a linear regression:

```r
mod_lm <- lm(LZINC ~ DAY, data = OYSTERZN)
sum_lm <- summary(mod_lm)
sum_lm
```

```
## 
## Call:
## lm(formula = LZINC ~ DAY, data = OYSTERZN)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.8974 -0.2448 -0.0709  0.2958  0.5715 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  6.073833   0.056834   106.9   <2e-16 ***
## DAY         -0.005314   0.000243   -21.9   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.354 on 71 degrees of freedom
## Multiple R-squared: 0.871,	Adjusted R-squared: 0.869 
## F-statistic:  478 on 1 and 71 DF,  p-value: <2e-16
```


which is fitting the model
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=ln%28Zn%29%20=%20a%20*%20day%2Bintercept) with a = `-0.0053` and intercept = `6.07`

Now plot data and model, as well as the residuals:

```r
# fitted values
fit_mod_lm <- fitted(mod_lm)

# data + fitted
plot(LZINC ~ DAY, OYSTERZN)
lines(OYSTERZN$DAY, fit_mod_lm)
```

![plot of chunk p85_lm_fitted](figure/p85_lm_fitted1.png) 

```r

# residuals
plot(mod_lm, which = 1)
```

![plot of chunk p85_lm_fitted](figure/p85_lm_fitted2.png) 


The mean square error can be calculated from the summary:

```r
# MSE
sum_lm$sigma^2
```

```
## [1] 0.1252
```

From which we can get an unbiased estimate of $C_0$:

```r
# unbiased estimate for C_0: exp(MSE/2) * exp(Intercept)
exp(sum_lm$sigma^2/2) * exp(sum_lm$coefficients[1, 1])
```

```
## [1] 462.4
```

where 

```r
sum_lm$coefficients[1, 1]
```

extracts the intercept from the summary.

The estimated k in the summary output is `-0.0053` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `2.4298 &times; 10<sup>-4</sup>`, and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_e%20=%20k%20-%20decayrate%20=%200.00531%20-%200.00283%20=%200.00248).

This result is similar to the weighted and non weighted nonlinear regression.
Again we have the same results as with SAS :) [Small deviations may be due to rounding error]




Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p85'.
### Quantitative Ecotoxicology, page 85, example 3.3, Backstripping

Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p89.csv) and read it into R:


```r
MERCURY <- read.table("p89.csv", header = TRUE, sep = ";")
```




```r
head(MERCURY)
```

```
##   DAY   PCI
## 1   3 27400
## 2   6 17601
## 3  10 12950
## 4  14 11157
## 5  22  9410
## 6  31  9150
```


```r
plot(PCI ~ DAY, data = MERCURY)
abline(v = 22)
```

![plot of chunk p89_raw](figure/p89_raw.png) 


We can identify two compartments on this plot: A slow one after day 22 and a fast one before day 22.

First we estimate ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_B) and ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_B) for the slow compartment, using linear regression of ln-transformed activity against day and predict from this slow compartment the activity over the whole period:


```r
# ln-transformation
MERCURY$LPCI <- log(MERCURY$PCI)
# fit linear model for day 31 to 94
mod_slow <- lm(LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
sum_mod_slow <- summary(mod_slow)
sum_mod_slow
```

```
## 
## Call:
## lm(formula = LPCI ~ DAY, data = MERCURY[MERCURY$DAY > 22, ])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2502 -0.0657  0.0605  0.0822  0.1386 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  9.43096    0.13987   67.42  2.6e-12 ***
## DAY         -0.01240    0.00213   -5.82   0.0004 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.135 on 8 degrees of freedom
## Multiple R-squared: 0.809,	Adjusted R-squared: 0.785 
## F-statistic: 33.9 on 1 and 8 DF,  p-value: 0.000395
```

```r
exp(coef(mod_slow)[1])
```

```
## (Intercept) 
##       12468
```

So this gives us the model ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C = 12468 e^{-0.0124 * Day}) for the slow component.

We do a bias correction and predict the activity for the whole data-range:

```r
# bias-correction
corr_mod_slow <- exp(sum_mod_slow$sigma^2/2)
# add bias corrected predicts to data.frame predict takes the whole
# data.frame es newdata, so we get predicts for every day.
MERCURY$PRED <- exp(predict(mod_slow, newdata = MERCURY)) * corr_mod_slow
# save C_B and k_B as objects (used later...)
CB <- exp(coef(mod_slow)[1]) * corr_mod_slow
KCB <- abs(coef(mod_slow)[2])
```


The residuals from these predictions for day 3 to 22 are associated with the fast compartment.
And we fit a linear regression to the ln-transformed residuals for the fast compartment.


```r
plot(LPCI ~ DAY, data = MERCURY)
points(MERCURY$DAY[1:4], MERCURY$LPCI[1:4], pch = 16)
abline(a = log(CB), b = -KCB)
for (i in 1:4) {
    lines(c(MERCURY$DAY[i], MERCURY$DAY[i]), c(log(MERCURY$PRED[i]), MERCURY$LPCI[i]), 
        lwd = 2)
}
```

![plot of chunk p89_residuals](figure/p89_residuals.png) 



```r
# extract residuals
MERCURY$ERROR <- MERCURY$PCI - MERCURY$PRED
MERCURY$ERROR[1:4]
```

```
## [1] 15276.20  5920.03  1834.41   579.42
```

```r
# fit linear model to ln(residuals) for Day 3 to 22
mod_fast <- lm(log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY < 22, ])
sum_mod_fast <- summary(mod_fast)
sum_mod_fast
```

```
## 
## Call:
## lm(formula = log(ERROR) ~ DAY, data = MERCURY[MERCURY$DAY < 22, 
##     ])
## 
## Residuals:
##       1       2       3       4 
##  0.0278 -0.0304 -0.0157  0.0183 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 10.49602    0.03755   279.5 0.000013 ***
## DAY         -0.29659    0.00407   -72.9  0.00019 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.0337 on 2 degrees of freedom
## Multiple R-squared:    1,	Adjusted R-squared: 0.999 
## F-statistic: 5.32e+03 on 1 and 2 DF,  p-value: 0.000188
```

```r
exp(coef(mod_fast)[1])
```

```
## (Intercept) 
##       36171
```


So the model for the fast component is: ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C = 36171 e^{-0.297 * Day})


```r
# bias correction
corr_mod_fast <- exp(sum_mod_fast$sigma^2/2)
# save C_A and k_A as objects
CA <- exp(coef(mod_fast)[1]) * corr_mod_fast
KCA <- abs(coef(mod_fast)[2])
```



Now we have two models: one for the fast component and one for the slow component, and we can make a plot similar to Figure 8.1 in Newman and Clements (2008, pp. 119–120). 


```r
plot(LPCI ~ DAY, data = MERCURY)
abline(mod_slow)
abline(mod_fast, lty = "dotted")
legend("topright", c("slow", "backstripped-fast"), lty = c("solid", "dashed"), 
    cex = 0.8)
```

![plot of chunk p89_backstripped](figure/p89_backstripped.png) 

```r
# Estimates
c(CA, KCA, CB, KCB)
```

```
## (Intercept)         DAY (Intercept)         DAY 
##  3.6192e+04  2.9659e-01  1.2583e+04  1.2403e-02
```


We can use this estimates as start-values to fit a non-linear Model to the data (therefore we stored them into objects).

We want to fit the following model:
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C=C_A * e^{-k_A*Day}%2BC_B e^{-k_B*Day})


```r
nls_mod1 <- nls(PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY), 
                data = MERCURY, 
                algorithm = "port",    # to use the bonds
                start = list(KCB = KCB, KCA = KCA, CB = CB, CA = CA),
                lower = c(0, 0, 5000, 20000), 
                upper = c(1, 1, 20000, 45000))
sum_nls_mod1 <- summary(nls_mod1)
sum_nls_mod1
```

```
## 
## Formula: PCI ~ CA * exp(-KCA * DAY) + CB * exp(-KCB * DAY)
## 
## Parameters:
##     Estimate Std. Error t value Pr(>|t|)    
## KCB 1.19e-02   1.41e-03    8.45  3.9e-06 ***
## KCA 2.97e-01   4.47e-02    6.66  3.6e-05 ***
## CB  1.22e+04   8.58e+02   14.27  1.9e-08 ***
## CA  3.79e+04   4.82e+03    7.86  7.7e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 701 on 11 degrees of freedom
## 
## Algorithm "port", convergence message: relative convergence (4)
```

* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_A) is estimated as `37909` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `4824`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=C_B) is estimated as `12245` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `858`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_A) is estimated as `0.297` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `0.045`
* ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=k_A) is estimated as `0.012` ![alt text](http://chart.apis.google.com/chart?cht=tx&chl=\\pm) `0.001`


And finally we plot data and model.

```r
plot(PCI ~ DAY, data = MERCURY, type = "n")
points(MERCURY$DAY, MERCURY$PCI, pch = ifelse(MERCURY$DAY <= 22, 16, 17))
# smooth line
pred_nls_mod1 <- predict(nls_mod1, newdata = data.frame(DAY = seq(0, 100, 1)))
lines(seq(0, 100, 1), pred_nls_mod1)
legend("topright", c("Fast", "slow"), pch = c(16, 17))
```

![plot of chunk p89_nls](figure/p89_nls.png) 



Again we get nearly the same results with R, except for some differences in the linear models.

This is probably due to the bias-correction in slow-component-model.
We have a MSE of

```r
sum_mod_slow$sigma^2
```

```
## [1] 0.018349
```

which is identical to the book. From the previous example, the bias can be estimated as
![alt text](http://chart.apis.google.com/chart?cht=tx&chl=e^{MSE/2}):

```r
exp(sum_mod_slow$sigma^2/2)
```

```
## [1] 1.0092
```

which is different to the book (1.002).


I have no SAS at hand, so I cannot check this with SAS. However let me know if there is an error
in my calculations.


**Refs**

> Newman, Michael C., and William Henry Clements. Ecotoxicology: A Comprehensive Treatment. Boca Raton: Taylor /& Francis, 2008. 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under file name 'p89'.
### Quantitative Ecotoxicology, Page 94, Example 3.5


Get the data from [here](https://raw.github.com/EDiLD/r-ed/master/quantitative_ecotoxicology/data/p94.csv) and read it into R:


```r
LEAD <- read.table("p94.csv", header = TRUE, sep = ";")
```




```r
head(LEAD)
```

```
##    DAY LEAD
## 1 0.16 41.0
## 2 0.16 31.0
## 3 0.16 25.3
## 4 1.00 30.5
## 5 1.00 22.7
## 6 1.00 22.0
```



As always we first take a look at the data:

```r
plot(LEAD ~ DAY, LEAD)
```

![plot of chunk p94_raw](figure/p94_raw.png) 


A simple power model may fit the data:

$$C_t = C_1~t^{−P}$$

We could fit such model as in example 3.3 via Nonlinear Least Squares or we could try to linearize the relationship by a ln-transform  of both DAY and LEAD:


```r
LEAD$LLEAD <- log(LEAD$LEAD)
LEAD$LDAY <- log(LEAD$DAY)
plot(LLEAD ~ LDAY, LEAD)
```

![plot of chunk p94_linear](figure/p94_linear.png) 


Now we can us lm() to estimate the coefficients and check our model:


```r
# fit model
mod <- lm(LLEAD ~ LDAY, data = LEAD)
```


The residuals show no pattern:

```r
plot(mod, which = 1)
```

![plot of chunk p94_residuals](figure/p94_residuals.png) 


From the model-output:

```r
mod_sum <- summary(mod)
mod_sum
```

```
## 
## Call:
## lm(formula = LLEAD ~ LDAY, data = LEAD)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.4568 -0.1789  0.0372  0.1689  0.4169 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   3.0008     0.0641   46.80  < 2e-16 ***
## LDAY         -0.2715     0.0313   -8.67  1.5e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.238 on 22 degrees of freedom
## Multiple R-squared: 0.773,	Adjusted R-squared: 0.763 
## F-statistic: 75.1 on 1 and 22 DF,  p-value: 1.53e-08
```



We see that out fitted model hast the formula:
$$Ln(LEAD) = 3.0008 - 0.272 ln(DAY)$$
with an R-squared of 0.77 and is statistically significant. The standard errors for the two parameters are 0.064 and 0.031.

So our backtransformed model would be:
$$ LEAD = exp(3.0008)~Day^{-0.272} = 20.68~Day^{-0.272}$$

Finally we can also plot our model:

```r
plot(LLEAD ~ LDAY, LEAD)
abline(mod)
```

![plot of chunk p94_model](figure/p94_model.png) 



Code and data are available at my [github-repo](https://github.com/EDiLD/r-ed/tree/master/quantitative_ecotoxicology) under filename 'p94'.


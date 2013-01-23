Today I added an interface to the **[freshwaterecology database](http://www.freshwaterecology.info)** to taxize.

<a href="http://www.freshwaterecology.info"><img src="http://www.freshwaterecology.info/images/freshwater_all.gif" alt="freshwater-header" width="400"/></a>

Currently only the macro-invertebrate database is supported, but I'll work also on the other four databases.

The functions have the ```fresh``` prefix and to use them you need the development version from github:


```r
install.packages("devtools")
require(devtools)
install_github("taxize_", "ropensci")
require(taxize)
```





Example
-------------------------

Let's say we have a species list:

```r
spec <- c("Aeshna grandis", "Baetis rodani", "Anax imperator", "Asellus aquaticus", 
    "Gammarus pulex")
```


We can check if our taxa are valid using the [Taxa Validation Tool (TVT)](http://www.freshwaterecology.info/TaxaDB_TVT.php). 

From within R we can use the new function ```fresh_validate```:

```r
spec_val <- fresh_validate(spec)
spec_val
```

```
##      Status    Genus   Species         Submitted
## 2     valid   Aeshna   grandis    Aeshna grandis
## 3 not valid                        Baetis rodani
## 4     valid     Anax imperator    Anax imperator
## 5     valid  Asellus aquaticus Asellus aquaticus
## 6     valid Gammarus     pulex    Gammarus pulex
```


We see that all taxa except the second (obviously a typo) are valid. 

We can also query the AQEM, DV, TCM, Furse, Perla, Ecoprof-Codes for these species (using ```fresh_codes```) or query ecological parameters (using ```fresh_traits``).

Currently the freshwaterecology-database comprises 40 Traits with 224 Modalities. You can get an overview from the description-table ```fresh_desc``` which is part of the taxize package:


```r
tail(fresh_desc)
```

```
##     Modality                    Description                          Trait
## 219      tec           terrestrial clutches                   reproduction
## 220      ase                        asexual                   reproduction
## 221      pas                      parasitic                   reproduction
## 222      rst                     r-strategy                  r- K-strategy
## 223      kst                     K-strategy                  r- K-strategy
## 224      oil occurrence in large quantities occurrence in large quantities
```


Getting the traits for our species is simple: We can use the object returned by ```fresh_validate``` (which contains a cookie) with ```fresh_traits```:


```r
spec_traits <- fresh_traits(spec_val)
```






```fresh_traits``` queries always all available ecological parameters, which yields to a very wide data.frame:

```r
dim(spec_traits)
```

```
## [1]   5 227
```


I we want to have only a subset of traits, we can use the above mentioned description table:


```r
take <- fresh_desc$Modality[fresh_desc$Trait %in% c("reproduction", "respiration")]
spec_traits[, c("x", take)]
```

```
##                   x teg gil pls spi ves tap sur ovo fie cie fic frc vec
## 1    Aeshna grandis   0   0   0   0   0   0   0   0   0   0   0   0   0
## 2     Baetis rodani  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
## 3    Anax imperator   0   0   0   0   0   0   0   0   0   0   0   0   0
## 4 Asellus aquaticus   0   0   0   0   0   0   0   0   0   0   0   0   0
## 5    Gammarus pulex   0   0   0   0   0   0   0   0   0   0   0   0   0
##   tec ase pas
## 1   0   0   0
## 2  NA  NA  NA
## 3   0   0   0
## 4   0   0   0
## 5   0   0   0
```


Integrating the other four databases is on my ToDo-list.


The fresh_* functions are under development, so if you have suggestions or find bugs contact me:
either on [github](https://github.com/ropensci/taxize_/issues?state=open) or [on this site](http://edild.github.com/contact/).

<a href="http://ropensci.org"><img src="http://ropensci.org/alm/images/ropensci.png" alt="ropensci" width="200"/></a>





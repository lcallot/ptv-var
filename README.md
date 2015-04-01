# Replication material for: 'Vector Autoregressions with Parsimoniously Time-Varying Parameters and an Application to Monetary Policy'.
## Laurent Callot and Johannes Tang Kristensen.


[Link to the paper](http://lcallot.github.io/papers/ptv-var/)

---

Date 01/04/2015

This repository contains the material necessary to replicate the __simulations__ results, the out from the empirical __application__, and the __admissible region__ figure in: 'Vector Autoregressions with Parsimoniously Time-Varying Parameters and an Application to Monetary Policy'. 

### Required packages 

The __parsimonious__ package from the eponymous repository is required:

```r
library('devtools')
install_github('lcallot/parsimonious')
library('parsimonious')
```

In addition the following packages should be installed : _ggplot2_, _gridExtra_, _sandwich_, _quantmod_, _xts_, _parallel_, and _xtable_, all available from CRAN.   


### /application

+ __PTVApp.R__ estimates the model with two lags and produces the figures _out.pdf_ and _inf.pdf_. 
+ __PTVApp\_1lag.R__ estimates the two models with 1 lag and produces _1lag.pdf_ and _ptv\_rho.pdf_.
+ __FREDData.RData__ contains the data downloaded from the St. Louis Fed. If _FREDData.RData_ is not found when running one of the two previous file the data is downloaded using FRED's API and automatically formatted. Use it if you want to automatically update the plots with more recent data.   
+ __Libs.R__ contains functions.


### /simulations

The _simulations_ folder contains 4 sub-folders. 

+ __scripts__ contains 4 R scripts to replicate the 4 experiments presented in the paper. The script run Monte Carlo iterations in parallel using _mclapply_. The scripts are configures to perform 10000 replications on 16 cores, make sure to adjust these settings. A file containing the simulation output is saved in _../mcsaves/_.
+ __mcsaves__ an empty folder, the simulation script saves R objects with the simulation results. Those are not included because of their size. 
+ __knitr__ contains a series of _.Rnw_ files to format the simulation output from _mcsaves_ and generate the plots and tables. These _knitr_ files are not meant to generate pretty PDFs, just useful _tex_ code for the tables and save the figures in pdf. 
+ __subs__ contains functions. 

### /admissible
The file __adm\_region.R__ folder contains the code to generate Figure 1 in the paper using _ggplot2_.  



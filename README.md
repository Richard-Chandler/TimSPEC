# TimSPEC: R package for postprocessing ensembles of climate time series projections

This package contains routines for the analysis of ensembles of time series of climate projections 
from multiple climate models, accounting for complex shared discrepancies between properties of the models
and the real climate system. Routines are provided for situations where the climate models can be considered 
as exchangeable, and also for structured ensembles using coupled pairs of models (such as global and regional 
model pairs in regional modelling experiments). The methodology is suitable for use with time series 
containing a single value per year. It is based around Gaussian dynamic linear models, although options are 
provided for ensuring that postprocessed projections remain non-negative if required. 

The methodology is described in detail in 

> Chandler, R.E., C.R. Barnes and C.M. Brierley (2023): _Decision-relevant characterisation of uncertainty in UK climate projections._ 
[Technical report, UK Climate Resilience Programme project CR20-3 _Enabling the use and producing improved understanding of EuroCORDEX data over the UK_](https://www.ucl.ac.uk/statistics/sites/statistics/files/projectionuncertainty.pdf).

# Requirements

The package has been tested under [R (version 4.1.2 and later)](https://www.r-project.org/), under both `Windows` and `Ubuntu` operating systems. 
It makes use of the `dlm`, `Hmisc`, `magick` and `numDeriv` add-on packages in `R`. **Note, however:** there are memory leaks in the `dlm` package as 
distributed via [CRAN](https://cloud.r-project.org/web/packages/dlm/index.html). The package author has been notified but has not fixed them: 
a patched version is needed, therefore. See the 'Installation' section below. 

# Installation

To install the package from the `R` command line:

```{r}
bibble
```

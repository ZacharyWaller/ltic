# Methods for Left Truncated Interval-Censored Survival (LTIC) data 

This package implements a fast product-limit algorithm to calculate the 
non-parametric maximum likelihood estimator for LTIC data. This also includes 
simulations and plots from the paper "A fast and stable NPMLE estimator for 
left-truncated and interval-censored data".

# Installation

```{r}
# install.packages("devtools")
devtools::install_github("ZacharyWaller/ltic")
```

# MHCPS Data

The above paper applies these methods to the Massachusetts Health Care Panel 
dataset. This data is available from the supporting information from Pan and 
Chappell (2002), [here](https://doi.org/10.1111/j.0006-341X.2002.00064.x).

Some changes need to be made to that data to make it make sense - there are some 
lines overlapping and lack of spaces between certain columns. Aside from the new 
lines changes, the other changes are shown in `tests/mhcps.R`.

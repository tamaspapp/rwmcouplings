This repository contains the R package `rwmcouplings`, which reproduces the numerical experiments in "A new and asymptotically optimally contracting coupling for the random walk Metropolis" (2022) by Tamas Papp and Chris Sherlock.

The package can be installed via `devtools`:
```
install.packages("devtools")
devtools::install_github("tamaspapp/rwmcouplings")
```

The `/inst/` directory contains all of the scripts to reproduce the experiments.
1. The scripts `spherical.R`, `elliptical.R` and `svm.R` reproduce the experiments in the main text.
2. The scripts `asymptoticallyoptimal.R` and `crnhugvis3d.R` reproduce additional visualisations in the appendices.

The scripts require additional external packages, which can be installed as follows:
```
# Data processing
install.packages("deSolve")
install.packages("mvtnorm")
install.packages("doParallel")
install.packages("doRNG")
install.packages("reshape2")

# Plotting
install.packages("ggplot2")
install.packages("ggh4x")
install.packages("latex2exp")
install.packages("ggpubr")
install.packages("gridExtra")

# Hug&Hop coupling visualisation
install.packages("geosphere")
```


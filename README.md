This repository contains the R package `rwmcouplings` that reproduces the numerical experiments of the paper "[Scalable couplings for the random walk Metropolis algorithm](https://arxiv.org/abs/2211.12585)" by Tamas P. Papp and Chris Sherlock.

The package can be installed via `devtools`:
```
install.packages("devtools")
devtools::install_github("tamaspapp/rwmcouplings")
```

The `inst/` directory contains all of the scripts to reproduce the experiments. For the main text:
1. `spherical.R` reproduces Figure 1.
2. `scaling/stepsizescaling-gcrn.R` reproduces Figure 2.
3. `scaling/asymptote.R` reproduces Figure 3.
4. `elliptical.R` reproduces Figure 4.
5. The scripts in `svm/` reproduce the experiments of Sections 7.1-7.3.
6. The scripts in `binreg/` reproduce the experiments of Section 7.4. Note that `binreg/binreg_coupling.R` was run on a server CPU with 56 cores.

Additional external packages may be required by the scripts, including:
```
install.packages("deSolve")
install.packages("doParallel")
install.packages("doRNG")
install.packages("mvtnorm")
install.packages("reshape2")

# Plotting
install.packages("ggh4x")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("gridExtra")
install.packages("latex2exp")

# Hug & Hop coupling visualisation from the appendix
install.packages("geosphere")
```

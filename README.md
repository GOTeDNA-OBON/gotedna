<!-- README.md is generated from README.Rmd. Please edit that file -->

# GOTeDNA

## An R package for guidance on optimal eDNA sampling periods to develop, optimize, and interpret monitoring programs

### This README is intended for installation and usage of the app. There is also a starter guide to the codebase provided [here](README_for_coders.md) in README_for_coders.md
<!-- badges: start -->

<!-- badges: end -->

The goal of GOTeDNA is to import and format eDNA qPCR and metabarcoding metadata/data from GOTeDNA sample templates, visualize species detection periods, and statistically delineate optimal species detection windows.

## Installation

### For non-R users

#### Install R

We recommend to use R and RStudio: <https://posit.co/download/rstudio-desktop/>

1.  Download R for your OS: <https://cran.rstudio.com/>

2.  Install R Studio

### Install the GOTeDNA package

#### Install Rtools (For Windows)

Note: macOS and Linux do not need Rtools to be installed to run the GOTeDNA package

The Rtools version appropriate for your R Version will need to be installed from https://cran.r-project.org/bin/windows/Rtools/ 

To see what R Version you currently have:

  ``` r
R.version.string
```

#### R users with access to the GitHub repository

##### Clone the repository with synced package versions (assumes git is already installed on your machine)

`git clone https://github.com/AnaisLacoursiereRoussel/GOTeDNA.git`

`cd GOTeDNA`

Start R in this folder, then in R, install renv if needed:

`install.packages("renv")`

Restore the project environment from the lockfile

`renv::restore()`

Follow terminal instructions for installing/syncing packages.

*Download and installation could take 10-20 minutes.

## Usage

To load the package: 

``` r
library("GOTeDNA")
```

### Shiny

The Shiny application can be launched with:

``` r
run_gotedna_app()
```

### Import data

To import your data within GOTeDNA, it must be formatted within the GOTeDNA template Excel sheets.  

Please refer to Appendices 1 and 2 in the GOTeDNA manuscript to access the sample metadata and metabarcoding templates.

### Visualization

The GOTeDNA app displays the following visualizations for each selected taxon and set of parameters (e.g. detection threshold, primer, etc.) 

#### Species monthly detection

![](man/figures/README-smooth-3.png)

#### Monthly detection probabilities
###### Predicted sample size required to reach targeted probability of detection given the month of sampling.
![](man/figures/README-thresh_fig-3.png)

#### Heat map
###### Heatmap displaying variation of normalized species eDNA detection probability.
![](man/figures/README-hm-3.png)

#### Effort needed

![](man/figures/README-effort-3.png)

#### Sampling effort
###### Monthly water sampling effort and proportion of samples having positive eDNA detection.

![](man/figures/README-field-4.png)

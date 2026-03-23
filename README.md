<!-- README.md is generated from README.Rmd. Please edit that file -->

# GOTeDNA

## An R package for guidance on optimal eDNA sampling periods to develop, optimize, and interpret monitoring programs

The goal of GOTeDNA is to import and format eDNA metabarcoding metadata/data from GOTeDNA sample templates, visualize species detection periods, and statistically delineate optimal species detection windows.

#### This README is intended for installation and usage of the app. There is also a starter guide to the codebase provided [here](README_for_coders.md) in README_for_coders.md
<!-- badges: start -->

<!-- badges: end -->

## For New-Coders
### Install R

Installation of R and RStudio are currently required to access the app.

This link will guide installation of both applications based on your operating system (Linux, Mac, Windows): <https://posit.co/download/rstudio-desktop/>

## Requirements for loading the GOTeDNA Application
### The following hints are provided to aid loading of the app. 

New coders will need to follow the hints specific to their device before installation of GOTeDNA. 

Users with previous coding experience may have these packages already loaded and can install the app directly by following instructions in section "Installation of the GOTeDNA Application".
If errors occur, refer to the hints below. 

#### Hint (for Windows)

The Rtools version appropriate for your R Version will need to be installed. 

To see what R version you currently have, paste the following code into the R Console:

``` r
R.version.string
```
Then follow the instructions from https://cran.r-project.org/bin/windows/Rtools/ 

#### Hint (for MacOS)

The package manager Homebrew will need to be installed. Instructions for installation can be found here: https://brew.sh/

Once installed, run the following line of code in the terminal:

``` r
brew install pkg-config proj geos gdal
```


#### Hint (For Linux)
Missing development libraries and compilers may need to be installed for R and RStudio. 

To help mitigate potential errors, read the warnings() that occur during the installation of pak and GOTeDNA to determine which libraries are missing.

#### Example for Debian/Devuan

Installing R
``` r
sudo apt install r-base r-base-dev
```

Installing RStudio
``` r
wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2026.01.1-403-amd64.deb
sudo apt install gdebi ## gdebi will fetch dependencies for rstudio if they are missing
sudo gdebi rstudio-2026.01.1-403-amd64.deb
```
 
Example of missing libraries (partial list)
``` r
sudo apt install libproj-dev proj-bin
sudo apt install libgdal-dev gdal-bin
```
Note: The specific command and package names differ among Linux Distributions

## Installation of the GOTeDNA Application
Install required packages by pasting the following code into the R Console:

``` r
install.packages("pak")

pak::pak(c(
  "trafficonese/leaflet.extras",
  "AnaisLacoursiereRoussel/GOTeDNA"
))
```

Load the GOTeDNA library: 

``` r
library(GOTeDNA)
```

Launch the Shiny application in a new browser: 

``` r
run_gotedna_app()
```
## Functionalities within the App
### Import data (Optional)

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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# GOTeDNA

## An R package for guidance on optimal eDNA sampling periods to develop, optimize, and interpret monitoring programs

<!-- badges: start -->

<!-- badges: end -->

The goal of GOTeDNA is to import and format eDNA qPCR and metabarcoding metadata/data from GOTeDNA sample templates, visualize species detection periods, and statistically delineate optimal species detection windows.

## Installation

### non-R users

#### Install R

We recommend to use R and RStudio: <https://posit.co/download/rstudio-desktop/>

1.  Download R for your OS: <https://cran.rstudio.com/>

2.  Install R Studio

#### Install the package

You first need to have access to the archive `GoteDNA_{version}.tar.gz`. Once you have obtained the archive, use the following:

``` r
install.packages("path/to/GOTeDNA_{version}.tar.gz")
```

### R users with access to the GitHub repository

You can install the development version of GOTeDNA from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AnaisLacoursiereRoussel/GOTeDNA", dependencies = TRUE)
```

Or if you have a local copy of the repo:

``` r
install.packages("devtools")
devtools::install_local("path/to/the/repo", dependencies = TRUE)
```

## Usage

### R function categories:

-   Import data
-   Clean/tidy data
-   Visualization
-   Shiny

``` r
library("GOTeDNA")
```

### Import data

To import your data within GOTeDNA, it must be formatted within the GOTeDNA template Excel sheets prior to calling in the `read_data()` function.

Please contact [Anais.Lacoursiere\@dfo-mpo.gc.ca](mailto:Anais.Lacoursiere@dfo-mpo.gc.ca) for access to the latest templates.

```         
D_mb_ex <- read_data(choose.method = "metabarcoding", path.folder = NULL)
```

### Clean/tidy data

As the Shiny app controls the filtering of each function internally, please filter to the taxonomy level and name that you wish to explore with `dplyr::filter()` when working with the code outside the app.

The example data herein contains a sample of metabarcoding data from a single protocol, and contains only the species detected within genus *Acartia*.

``` r
newprob <- calc_det_prob(
  data = D_mb_ex |> dplyr::filter(genus == "Acartia")
  )

scaledprobs <- scale_newprob(
  data = D_mb_ex |> dplyr::filter(genus == "Acartia"), 
  newprob
  )

win <- calc_window(
  threshold = "75",
  scaledprobs = scaledprobs |> dplyr::filter(species == "Acartia longiremis")
  )

win$opt_sampling
win$fshTest
```

### Visualization

#### Species monthly detection

``` r
 smooth_fig(
   data = D_mb_ex |> dplyr::filter(species == "Acartia longiremis")
   )
```

![](man/figures/README-smooth-1.png)

#### Monthly detection probabilities

``` r
thresh_fig(
  threshold = "75",   
  scaledprobs = scaledprobs |> dplyr::filter(species == "Acartia longiremis")
)
```

![](man/figures/README-thresh_fig-1.png)

#### Heat map

``` r
hm_fig(
  scaledprobs = scaledprobs |> dplyr::filter(class == "Copepoda")
)
```

![](man/figures/README-hm-1.png)

#### Effort needed

``` r
effort_needed_fig(
   scaledprobs = scaledprobs |> dplyr::filter(species == c("Acartia longiremis","Acartia tonsa"))
)
```

![](man/figures/README-effort-1.png)

#### Sampling effort

``` r
field_sample_fig(
  data = D_mb_ex |> dplyr::filter(class == "Copepoda")
)
```

![](man/figures/README-field-1.png)

### Shiny

Once the package is loaded, the Shiny application can be launched using the following function:

``` r
run_gotedna_app()
```

There is also a docker file available:

``` sh
# build the container 
docker build -t gotedna . 
# use the container
docker run -it --rm --network host gotedna 
# the shiny will be available at http://0.0.0.0:9292
```

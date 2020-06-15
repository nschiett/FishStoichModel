[![DOI](https://zenodo.org/badge/256277906.svg)](https://zenodo.org/badge/latestdoi/256277906)



# Nutrient limitation, bioenergetics, and stoichiometry: a new model to predict elemental fluxes mediated by fishes

This repository contains code and data needed to reproduce the figures and tables of the manuscript:

Schiettekatte, N. M. D., Barneche, D. R., Villéger, S., Allgeier, J, Burkepile, D.; Brandl, S.J., Casey, J. M.; Mercière, A; Munsterman, K; Morat, F; Parravicini, V (2020) (in press.) Nutrient limitation, bioenergetics, and stoichiometry: a new model to predict elemental fluxes mediated by fishes. *Functional ecology*

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies = TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies = TRUE)
```

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

To replicate the analyses, it is nescessary to install `fishflux` through GitHub. To do so, follow instructions in this [link](https://github.com/nschiett/fishflux).

Then, to generate all figures and tables, simply run:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that all the combined Bayesian models in this paper will take a several hours (up to a day) to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename = 'build.R')
```

### This paper was produced using the following software and associated packages:
```
R version 3.6.2 (2019-12-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] egg_0.4.5           gridExtra_2.3       tidybayes_2.0.1     knitr_1.28          rmarkdown_2.1       plyr_1.8.5          png_0.1-7           fishualize_0.2.0   
 [9] dplyr_0.8.5         purrr_0.3.3         ggplot2_3.3.0       fishflux_0.0.0.9001 Rcpp_1.0.4.6       

loaded via a namespace (and not attached):
 [1] lattice_0.20-40       tidyr_1.0.2           prettyunits_1.1.1     ps_1.3.2              assertthat_0.2.1      digest_0.6.25         R6_2.4.1             
 [8] stats4_3.6.2          coda_0.19-3           evaluate_0.14         httr_1.4.1            pillar_1.4.3          rlang_0.4.5           rstudioapi_0.11      
[15] callr_3.4.3           labeling_0.3          stringr_1.4.0         loo_2.2.0             munsell_0.5.0         xfun_0.12             compiler_3.6.2       
[22] rstan_2.19.3          pkgconfig_2.0.3       pkgbuild_1.0.6        htmltools_0.4.0       rstantools_2.0.0.9000 tidyselect_1.0.0      tibble_3.0.0         
[29] arrayhelpers_1.1-0    codetools_0.2-16      matrixStats_0.55.0    fansi_0.4.1           remake_0.3.0          crayon_1.3.4          withr_2.1.2          
[36] gtable_0.3.0          lifecycle_0.2.0       magrittr_1.5          storr_1.2.1           StanHeaders_2.19.2    scales_1.1.0          cli_2.0.2            
[43] stringi_1.4.6         farver_2.0.3          ellipsis_0.3.0        vctrs_0.2.4           cowplot_1.0.0         tools_3.6.2           forcats_0.4.0        
[50] svUnit_0.7-12         glue_1.4.0            processx_3.4.2        parallel_3.6.2        yaml_2.2.1            inline_0.3.15         colorspace_1.4-1         
```

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

## Bug reporting
* Please [report any issues or bugs](https://github.com/nschiett/FishStoichModel/issues).

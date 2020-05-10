
<!-- README.md is generated from README.Rmd. Please edit that file -->
brainflowprobes
===============

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable) [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/brainflowprobes.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/brainflowprobes) [![BioC dev status](http://www.bioconductor.org/shields/build/devel/bioc/brainflowprobes.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/brainflowprobes) [![Codecov test coverage](https://codecov.io/gh/LieberInstitute/brainflowprobes/branch/master/graph/badge.svg)](https://codecov.io/gh/LieberInstitute/brainflowprobes?branch=master) [![R build status](https://github.com/LieberInstitute/brainflowprobes/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/LieberInstitute/brainflowprobes/actions) [![Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.](http://www.bioconductor.org/shields/posts/brainflowprobes.svg)](https://support.bioconductor.org/t/brainflowprobes/) [![GitHub issues](https://img.shields.io/github/issues/LieberInstitute/brainflowprobes)](https://github.com/LieberInstitute/brainflowprobes/issues) <!-- badges: end -->

This package can be used to annotate candidate target probe sequences designed to label RNA in nuclei isolated from human postmortem brain for flow cytometry. This package also includes functions that visualize candidate sequence expression in four postmortem brain datasets that address nuclear expression, effects of RNA degradation and flow cytometric sorting, and cell type-specific expression to increase the likelihood of success as a custom target probe. Please refer to the BrainFlow publication or the [PrimeFlow RNA(TM) Assay Kit](https://www.thermofisher.com/order/catalog/product/88-18005-204) literature for more information.

Documentation
-------------

For more information about `brainflowprobes` check the vignettes [through Bioconductor](http://bioconductor.org/packages/brainflowprobes) or at the [documentation website](http://lieberinstitute.github.io/brainflowprobes).

Installation instructions
-------------------------

Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `brainflowprobes` using from [Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("brainflowprobes")
```

Citation
--------

Below is the citation output from using `citation('brainflowprobes')` in R. Please run this yourself to check for any updates on how to cite **brainflowprobes**.

``` r
print(citation('brainflowprobes'), bibtex = TRUE)
#> 
#> Price AJ (2020). _Plots and annotation for choosing BrainFlow target
#> probe sequence_. doi: 10.18129/B9.bioc.brainflowprobes (URL:
#> https://doi.org/10.18129/B9.bioc.brainflowprobes),
#> https://github.com/LieberInstitute/brainflowprobes - R package version
#> 1.1.2, <URL: http://www.bioconductor.org/packages/brainflowprobes>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Plots and annotation for choosing BrainFlow target probe sequence},
#>     author = {Amanda J Price},
#>     year = {2020},
#>     url = {http://www.bioconductor.org/packages/brainflowprobes},
#>     note = {https://github.com/LieberInstitute/brainflowprobes - R package version 1.1.2},
#>     doi = {10.18129/B9.bioc.brainflowprobes},
#>   }
```

Please note that the `brainflowprobes` was only made possible thanks to many other R and bioinformatics software authors, which are cited either in the vignettes and/or the paper(s) describing this package.

Code of Conduct
---------------

Please note that the derfinderPlot project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

Development tools
-----------------

-   Continuous code testing is possible thanks to [GitHub actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/) through *[usethis](https://CRAN.R-project.org/package=usethis)*, *[remotes](https://CRAN.R-project.org/package=remotes)*, *[sysreqs](https://github.com/r-hub/sysreqs)* and *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized to use [Bioconductor's docker containers](https://www.bioconductor.org/help/docker/) and *[BiocCheck](https://bioconductor.org/packages/3.11/BiocCheck)*.
-   Code coverage assessment is possible thanks to [codecov](https://codecov.io/gh) and *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://lieberinstitute.github.io/brainflowprobes) is automatically updated thanks to *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to *[devtools](https://CRAN.R-project.org/package=devtools)* and *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

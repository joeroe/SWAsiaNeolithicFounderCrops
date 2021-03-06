
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The “Neolithic Founder Crops” in Southwest Asia: Research Compendium

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joeroe/SWAsiaNeolithicFounderCrops/main?urlpath=rstudio)

This repository contains the data and code for our paper:

> Arranz-Otaegui, Amaia, and Roe, Joe (in prep). *Revisiting the concept
> of the “Neolithic Founder Crops” in southwest Asia*. To appear in
> Vegetation History & Archaeobotany. <https://doi.org/xxx/xxx>

<!--
Our pre-print is online here:

> Authors, (YYYY). _Revisiting the concept of the "Neolithic Founder Crops" in southwest Asia_. Name of journal/book, Accessed 20 Jan 2022. Online at <https://doi.org/xxx/xxx>
-->

## Contents

The `analysis` directory contains:

-   [:file_folder: SI1.Rmd](/analysis/SI1.Rmd): R Markdown source code
    for the supplement to our manuscript describing the quantitative
    analysis. It includes code to reproduce the figures and tables
    generated by the analysis.
    -   [:file_folder: SI1.html](/analysis/SI1.html): rendered version
        of the supplement in HTML format.
    -   [:file_folder: SI1.pdf](/analysis/SI1.pdf): rendered version of
        the supplement in PDF format.
-   [:file_folder: data](/analysis/data): Data used in the analysis.
    -   [:file_folder: raw_data](/analysis/data/raw_data): Data from
        other sources – see manuscript for attribution.
    -   [:file_folder: derived_data](/analysis/data/derived_data):
        Collated and processed dataset used for our analysis.
-   [:file_folder: figures](/analysis/figures): Individual figures and
    tables generated in the analysis.

Additional R functions supporting the analysis are contained in the `R`
directory.

## Run the compendium locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the compendium as a zip from from this URL:
[master.zip](/archive/master.zip). After unzipping: - open the `.Rproj`
file in RStudio - run `devtools::install()` to ensure you have the
packages this analysis depends on (also listed in the
[DESCRIPTION](/DESCRIPTION) file). - finally, open
`analysis/analysis.Rmd` and ‘knit’ or run
`rmarkdown::render("analysis/analysis.Rmd")` in the R console to produce
the `analysis.html` and `analysis.pdf`.

## Citation

In addition to our paper above, please cite this compendium as:

> Roe, Joe and Arranz-Otaegui, Amaia, (2022). *Compendium of R code and
> data for ‘Revisiting the concept of the “Neolithic Founder Crops” in
> southwest Asia’*. Accessed 20 Jan 2022. Available at
> <https://doi.org/xxx/xxx>

## Licenses

**Text and figures:**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code:** See the [DESCRIPTION](DESCRIPTION) file

**Data:** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

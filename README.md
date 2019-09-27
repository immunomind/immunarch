[![Follow](https://img.shields.io/twitter/follow/immunomind.svg?style=social)](https://twitter.com/intent/follow?screen_name=immunomind)
[![CI](https://gitlab.com/immunomind/immunarch/badges/master/pipeline.svg?style=flat-square)](https://gitlab.com/immunomind/immunarch/-/jobs)
[![Issues](https://img.shields.io/github/issues/immunomind/immunarch?style=flat-square)](http://github.com/immunomind/immunarch/issues)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3367200.svg)](https://doi.org/10.5281/zenodo.3367200)


# `immunarch` - An R Package for Painless Analysis of Large-Scale Immune Repertoire Data

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Installation Troubleshooting](#installation-troubleshooting)
- [Quick Start](#quick-start)
- [Bugs and Issues](#bugs-and-issues)
- [Citation](#citation)


## Introduction
`immunarch` is an R package designed to analyse TCR and BCR (immunoglobulin) repertoires, which constitute a large amount of data. The mission of `immunarch` is to make immune sequencing data analysis as effortless as possible and help you focus on research instead of coding. Follow us on [Twitter](https://twitter.com/immunomind) or [Telegram](https://t.me/immunomind) for news and updates.

## Features
1. Fast and easy manipulation of immune repertoire data:

    + The package automatically detects the format of your files---no more guessing what format is *that* file, just pass them to the package;
  
    + Supports all popular TCR and BCR analysis and post-analysis formats: [ImmunoSEQ](https://www.adaptivebiotech.com/products-services/immunoseq/), [IMGT](http://www.imgt.org/IMGTindex/IMGTHighV-QUEST.php), [MiTCR](https://github.com/milaboratory/mitcr/), [MiXCR](https://milaboratory.com/software/mixcr/), [MiGEC](https://milaboratory.com/software/migec/), [MigMap](https://github.com/mikessh/migmap), [VDJtools](https://milaboratory.com/software/vdjtools/), [tcR](https://github.com/imminfo/tcr), [AIRR](http://docs.airr-community.org/en/latest/), [10XGenomics](https://support.10xgenomics.com/single-cell-vdj/datasets/), [ArcherDX](https://archerdx.com/immunoverse). More coming in the future;

    + Works on any data source you are comfortable with: R data frames, data tables from [data.table](http://r-datatable.com), databases like [MonetDB](https://github.com/MonetDB), Apache Spark data frames via [sparklyr](https://spark.rstudio.com/);
    
    + Tutorial is available [here](https://immunarch.com/articles/v2_data.html).

2. Immune repertoire analysis made simple:

    + Most methods are incorporated in a couple of main functions with clear naming---no more remembering tens and tens of functions with obscure names. For details see [link](https://immunarch.com/articles/v3_basic_analysis.html);

    + Repertoire overlap analysis *(common indices including overlap coefficient, Jaccard index and Morisita's overlap index)*. Tutorial is available [here](https://immunarch.com/articles/v4_overlap.html);
  
    + Gene usage estimation *(correlation, Jensen-Shannon Divergence, clustering)*. Tutorial is available [here](https://immunarch.com/articles/v5_gene_usage.html);

    + Diversity evaluation *(ecological diversity index, Gini index, inverse Simpson index, rarefaction analysis)*. Tutorial is available [here](https://immunarch.com/articles/v6_diversity.html);

    + Tracking of clonotypes;
    
    + Coming in the next releases: CDR3 amino acid physical and chemical properties assessment, Kmer distribution measures and statistics, mutation networks.

3. Publication-ready plots with a built-in tool for visualisation manipulation: 

    + Rich visualisation procedures with [ggplot2](https://ggplot2.tidyverse.org/);
  
    + Built-in tool `FixVis` makes your plots publication-ready: easily change font sizes, text angles, titles, legends and many more with clear-cut GUI;
    
    + Tutorial is available [here](https://immunarch.com/articles/v7_fixvis.html).

## Installation
You can find the list of releases of immunarch here: https://github.com/immunomind/immunarch/releases

In order to install immunarch, you need to download it first. If you want to download the latest version, you need to download the package file, available here https://github.com/immunomind/immunarch/raw/master/immunarch.tar.gz

Note that you should not un-archive it!

After downloading the file, you need to install a number of packages with R commands listed below, and run the newly installed `devtools` package to install `immunarch` locally. Upon completion the dependencies will have been already downloaded and installed.
```r
install.packages("devtools", dependencies = T)
devtools::install_local("path/to/your/folder/with/immunarch.tar.gz", dependencies=T)
```

That's it, you can start using `immunarch` now!

## Installation troubleshooting
If you can not install dependencies, please try manual installation of all dependencies by executing the following command in R console.
```r
install.packages(c("rematch", "prettyunits", "forcats", "cellranger", "progress", "zip", "backports", "ellipsis", "zeallot", "SparseM", "MatrixModels", "sp", "haven", "curl", "readxl", "openxlsx", "minqa", "nloptr", "RcppEigen", "utf8", "vctrs", "carData", "pbkrtest", "quantreg", "maptools", "rio", "lme4", "labeling", "munsell", "cli", "fansi", "pillar", "viridis", "car", "ellipse", "flashClust", "leaps", "scatterplot3d", "modeltools", "DEoptimR", "digest", "gtable", "lazyeval", "rlang", "scales", "tibble", "viridisLite", "withr", "assertthat", "glue", "magrittr", "pkgconfig", "R6", "tidyselect", "BH", "plogr", "purrr", "ggsci", "cowplot", "ggsignif", "polynom", "fastcluster", "plyr", "abind", "dendextend", "FactoMineR", "mclust", "flexmix", "prabclus", "diptest", "robustbase", "kernlab", "GlobalOptions", "shape", "colorspace", "stringi", "hms", "clipr", "crayon", "httpuv", "mime", "jsonlite", "xtable", "htmltools", "sourcetools", "later", "promises", "gridBase", "RColorBrewer", "yaml", "ggplot2", "dplyr", "dtplyr", "dbplyr", "data.table", "gridExtra", "ggpubr", "heatmap3", "ggrepel", "reshape2", "DBI", "factoextra", "fpc", "circlize", "tidyr", "Rtsne", "readr", "shiny", "shinythemes", "treemap", "igraph", "airr", "ggseqlogo", "UpSetR", "stringr", "ggalluvial", "Rcpp"))
```


If you run in any other trouble, try the following steps:

1. Check your R version. Run `version` command in the console to get your R versions. If the R version is below 3.5.0 (for example, `R version 3.1.0`), try updating your R version to the latest one.

2. Check if your packages are outdated and update them. In RStudio you can run the "Update" button on top of the package list. In R console you can run the `old.packages()` command to view a list of outdated packages.

3. If you are on Mac and have issues like old packages can't be updated, or error messages such as `ld: warning: directory not found for option` or `ld: library not found for -lgfortran`, [this link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/) will help you to fix the issue.

4. If you are on Mac Mojave (1.14) and run into the following error:
```
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/math.h:301:15: fatal error: 'math.h' file not found
#include_next <math.h>
              ^~~~~~~~
```

Open Terminal, execute the following command and try again to install `immunarch`:
```
sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
```

5. If you are working under Linux and have issues with igraph library or have 
Fortran errors, see [this link](https://ashokragavendran.wordpress.com/2017/10/24/error-installing-rigraph-unable-to-load-shared-object-igraph-so-libgfortran-so-4-cannot-open-shared-object-file-no-such-file-or-directory/)

6. If you encounter the following error while running the `devtools::install_local` function:

```
1: In normalizePath(path.expand(path), winslash, mustWork) :
  path[1]="path/to/your/folder/with/immunarch.tar.gz":

In file.copy(x$path, bundle, recursive = TRUE) :
  problem copying No such file or directory
```

Check your path to the downloaded package archive file. It should not be "path/to/your/folder/with/immunarch.tar.gz", but a path on your PC to the downloaded file, e.g., "C:/Users/UserName/Downloads/immunarch.tar.gz" or ""/Users/UserName/Downloads/immunarch.tar.gz"".

7. If you are working under Windows and have issues with the package installation, or if you want to change the folder for R packages, feel free to check [this forum post](https://community.rstudio.com/t/help-regarding-package-installation-renviron-rprofile-r-libs-r-libs-site-and-r-libs-user-oh-my/13888/8).

8. If troubles still persist, message us on support@immunomind.io or create an issue in https://github.com/immunomind/immunarch/issues with the code that represents the issue and the output you get in the console.

## Quick Start
Importing data into R is fairly simple. The gist of the typical TCR or BCR explorational data analysis workflow can be reduced to the next few lines of code:
```r
# Load the package
library(immunarch)

# Load the data to the package
immdata = repLoad("path/to/your/folder/with/repertoires")
# If you folder contains metadata.txt file, immdata will have two elements:
# - immdata$data with a list of parsed repertoires
# - immdata$meta with the metadata file

# If you don't have your data, you can load repertoire data
# that comes with immunarch. Uncomment (i.e., remove the "#" symbol)
# in the next line to load the data
# data(immdata)

# Compute and visualise overlap statistics
ov = repOverlap(immdata$data)
vis(ov)

# Cluster samples using K-means algorithm applied to the number of overlapped clonotypes
# and visualise the results
ov.kmeans = repOverlapAnalysis(ov, .method = "mds")
vis(ov.kmeans)

# Compute and visualise gene usage with samples, grouped by their disease status
gu = geneUsage(immdata$data)
vis(gu, .by="Status", .meta=immdata$meta)

# Compute Jensen-Shannon divergence among gene distributions of samples, 
# cluster samples using the hierarchical clustering and visualise the results
gu.clust = geneUsageAnalysis(gu, .method = "js+hclust")
vis(gu.clust)

# Compare diversity of repertoires and visualise samples, grouped by two parameters
div = repDiversity(immdata$data, .method = "chao1")
vis(div, .by=c("Status", "Lane"), .meta=immdata$meta)

# Manipulate the visualisation of diversity estimates to make the plot publication-ready
div.plot = vis(div, .by=c("Status", "Lane"), .meta=immdata$meta)
fixVis(div.plot)
```

If you want to test the package without parsing any data, you can load a small test dataset provided along with the package. Load the data with the following command:

```r
data(immdata)
```
# Bugs and Issues

The mission of Immunarch is to make immune repertoires painless to analyse. All bug reports, documentation improvements, enhancements and ideas are welcome.

If through using immunarch you have an idea of your own or are looking for something in the documentation and thinking ‘this can be improved’...you can do something about it! Just let us know.

Bug reports are an important part of making immunarch more stable. Having a complete bug report will allow us to reproduce the bug and provide insight into fixing.

Bug reports must: 

1. Include a short, self-contained R snippet reproducing the problem. 
2. Add minimal data sample for us to reproduce the problem. If for some reasons you don't want to share it publicly on Gihub we are always available through [support@immunomind.io](mailto:support@immunomind.io).
3. Explain why the current behavior is wrong/not desired and what you expect instead.
4. If the issue is somehow connected with plotting or visualization, please attach a picture. It'll be much simple for us to see what you see.

# Citation

ImmunoMind Team. (2019). immunarch: An R Package for Painless Analysis of Large-Scale Immune Repertoire Data. Zenodo. http://doi.org/10.5281/zenodo.3367200

BibTex:
```
@misc{immunomind_team_2019_3367200,
  author       = {{ImmunoMind Team}},
  title        = {{immunarch: An R Package for Painless Analysis of 
                   Large-Scale Immune Repertoire Data}},
  month        = aug,
  year         = 2019,
  doi          = {10.5281/zenodo.3367200},
  url          = {https://doi.org/10.5281/zenodo.3367200}
}
```

For EndNote citation import the [`immunarch-citation.xml`](https://github.com/immunomind/immunarch/releases/download/latest/immunarch-citation.xml) file.

# License

The package is freely distributed under the Apache v2 license. You can read more about it [here](https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)).

Additionally, we provide an annual subscription that includes next services:

- 100 hours of consultations on the TCR & BCR repertoire analysis (contact us to purchase more);
- Priority email and call support;
- Package modifications and feature implementations are issued promptly; 
- Setup a cloud or cluster installation of *immunarch*, including the development of cloud *immunarch*-based software;
- Use *immunarch* team expertise in your projects;
- If you need license other than the current, contact us.

Contact us at support@immunomind.io for more information.

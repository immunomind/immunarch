[![Follow](https://img.shields.io/twitter/follow/immunomind.svg?style=social)](https://twitter.com/intent/follow?screen_name=immunomind)
[![CI](https://gitlab.com/immunomind/immunarch/badges/master/pipeline.svg?style=flat-square)](https://gitlab.com/immunomind/immunarch/-/jobs)
[![Issues](https://img.shields.io/github/issues/immunomind/immunarch?style=flat-square)](http://github.com/immunomind/immunarch/issues)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3367200.svg)](https://doi.org/10.5281/zenodo.3367200)


# `immunarch` - An R Package for Painless Bioinformatics Analysis of T-cell and B-cell Immune Repertoire Data

- [Introduction](#introduction)
- [Installation](#installation)
- [Features](#features)
- [Quick Start](#quick-start)
- [Bugs and Issues](#bugs-and-issues)
- [Citation](#citation)


## Introduction

`immunarch` is an R package designed to analyse T-cell receptor (TCR) and B-cell receptor (BCR) repertoires, aimed at medical scientists and bioinformaticians. The mission of `immunarch` is to make immune sequencing data analysis as effortless as possible and help you focus on research instead of coding. Follow us on [Twitter](https://twitter.com/immunomind) for news and updates.


## Installation
In order to install `immunarch` execute the following R commands:

```r
install.packages("devtools") # skip this if you already installed devtools
devtools::install_url("https://github.com/immunomind/immunarch/raw/master/immunarch.tar.gz")
```

That's it, you can start using `immunarch` now! See the [Quick Start](#quick-start) section below to dive into immune repertoire data analysis. If you run in any trouble with installation, take a look at the [Installation Troubleshooting](https://immunarch.com/articles/v1_introduction.html#installation-troubleshooting) section.

Note that there are quite a lot of dependencies to install with the package because it installs all the widely-used packages for data analysis and visualisation. You got both the AIRR data analysis framework and Data Science package eco-system with only one command!

<!--
### Pre-release version installation

Since releasing on CRAN is limited to one release per one-two months, you can install the latest pre-release version with bleeding edge features and optimisations directly from a code repository. In order to install the latest pre-release version, you need to execute only two commands:

```r
install.packages("devtools") # skip this if you already installed devtools
devtools::install_url("https://github.com/immunomind/immunarch/raw/master/immunarch.tar.gz")
```
-->

You can find the list of releases of `immunarch` here: https://github.com/immunomind/immunarch/releases


## Features

1. Fast and easy manipulation of immune repertoire data:

    + The package automatically detects the format of your files---no more guessing what format is *that* file, just pass them to the package;
  
    + Supports all popular TCR and BCR analysis and post-analysis formats, including single-cell data: [ImmunoSEQ](https://www.adaptivebiotech.com/products-services/immunoseq/), [IMGT](http://www.imgt.org/IMGTindex/IMGTHighV-QUEST.php), [MiTCR](https://github.com/milaboratory/mitcr/), [MiXCR](https://milaboratory.com/software/mixcr/), [MiGEC](https://milaboratory.com/software/migec/), [MigMap](https://github.com/mikessh/migmap), [VDJtools](https://milaboratory.com/software/vdjtools/), [tcR](https://github.com/imminfo/tcr), [AIRR](http://docs.airr-community.org/en/latest/), [10XGenomics](https://support.10xgenomics.com/single-cell-vdj/datasets/), [ArcherDX](https://archerdx.com/immunology/). More coming in the future;

    + Works on any data source you are comfortable with: R data frames, data tables from [data.table](http://r-datatable.com), databases like [MonetDB](https://github.com/MonetDB), Apache Spark data frames via [sparklyr](https://spark.rstudio.com/);
    
    + Tutorial is available [here](https://immunarch.com/articles/v2_data.html).

2. Immune repertoire analysis made simple:

    + Most methods are incorporated in a couple of main functions with clear naming---no more remembering tens and tens of functions with obscure names. For details see [link](https://immunarch.com/articles/v3_basic_analysis.html);

    + Repertoire overlap analysis *(common indices including overlap coefficient, Jaccard index and Morisita's overlap index)*. Tutorial is available [here](https://immunarch.com/articles/web_only/v4_overlap.html);
  
    + Gene usage estimation *(correlation, Jensen-Shannon Divergence, clustering)*. Tutorial is available [here](https://immunarch.com/articles/web_only/v5_gene_usage.html);

    + Diversity evaluation *(ecological diversity index, Gini index, inverse Simpson index, rarefaction analysis)*. Tutorial is available [here](https://immunarch.com/articles/web_only/v6_diversity.html);

    + Tracking of clonotypes across time points, widely used in vaccination and cancer immunology domains. Tutorial is available [here](https://immunarch.com/articles/web_only/v8_tracking.html);
    
    + Kmer distribution measures and statistics. Tutorial is available [here](https://immunarch.com/articles/web_only/v9_kmers.html);
    
    + Coming in the next releases: CDR3 amino acid physical and chemical properties assessment, mutation networks.

3. Publication-ready plots with a built-in tool for visualisation manipulation: 

    + Rich visualisation procedures with [ggplot2](https://ggplot2.tidyverse.org/);
  
    + Built-in tool `FixVis` makes your plots publication-ready: easily change font sizes, text angles, titles, legends and many more with clear-cut GUI;
    
    + Tutorial is available [here](https://immunarch.com/articles/web_only/v7_fixvis.html).
    
    
## Quick start
The gist of the typical TCR or BCR data analysis workflow can be reduced to the next few lines of code.

**1) Load the package and the data**

```r
# 1.1) Load the package into R:
library(immunarch)

# 1.2a) To quickly test immunarch, load the test dataset:
data(immdata)

# 1.2b) To try immunarch on your own data, use the `repLoad` function on your data folder:
immdata = repLoad("path/to/your/folder/with/repertoires")
```

**2) Analyse repertoire similarity at the clonotype level**

```r
# 2.1) Find the number of shared clonotypes and visualise it:
ov = repOverlap(immdata$data)
vis(ov)

# 2.2) Cluster samples by their similarity:
ov.kmeans = repOverlapAnalysis(ov, .method = "mds+kmeans")
vis(ov.kmeans)
```

**3) Find repertoire differences in the Variable gene usage**

```r
# 3.1) Compute V gene usage and and highlight gene differences in groups with different clinical status:
gu = geneUsage(immdata$data)
vis(gu, .by="Status", .meta=immdata$meta)

# 3.2) Cluster samples by their V gene usage similarity:
gu.clust = geneUsageAnalysis(gu, .method = "js+hclust")
vis(gu.clust)
```

**4) Find differences in the diversity of repertoires**

```r
# 4.1) Compare diversity of repertoires and visualise samples, grouped by both clinical status and sequencing Lane:
div = repDiversity(immdata$data, .method = "chao1")
vis(div, .by=c("Status", "Lane"), .meta=immdata$meta)
```

**5) Manipulate plots to make them publication-ready**

```r
# 5.1) Manipulate the visualisation of diversity estimates to make the plot publication-ready:
div = repDiversity(immdata$data, .method = "chao1")
div.plot = vis(div, .by=c("Status", "Lane"), .meta=immdata$meta)
fixVis(div.plot)
```

**6) Advanced methods**

For advanced methods such as clonotype tracking, kmer analysis and public repertoire analysis see "Tutorials".


# Bugs and Issues

The mission of `immunarch` is to make immune repertoires painless to analyse. All bug reports, documentation improvements, enhancements and ideas are welcome.

If through using `immunarch` you have an idea of your own or are looking for something in the documentation and thinking 'this can be improved'... you can do something about it! Just let us know via [GitHub](https://github.com/immunomind/immunarch/issues) or [support@immunomind.io](mailto:support@immunomind.io).

Bug reports are an important part of making `immunarch` more stable. Having a complete bug report will allow us to reproduce the bug and provide insight into fixing.

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

For EndNote citation import the [`immunarch-citation.xml`](https://gitlab.com/immunomind/immunarch/raw/master/immunarch-citation.xml?inline=false) file.

Preprint on BioArxiv is coming soon.


# License

The package is freely distributed under the Apache v2 license. You can read more about it [here](https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)).

If you are interested in our commercial-grade platform for biomarker discovery and AIRR analysis, contact us via support@immunomind.io

<!--
- Package modifications and feature implementations are issued promptly;
- Use *immunarch* team expertise in your projects;
- Priority email and call support;
- 100+ hours of consultations on the TCR & BCR repertoire analysis;
- Setup a cloud or cluster installation of *immunarch*, including the development of cloud *immunarch*-based software;
- If you need license other than the current, contact us.

Contact us at support@immunomind.io for more information.
-->

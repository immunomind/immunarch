[![Follow](https://img.shields.io/twitter/follow/immunomind.svg?style=social)](https://twitter.com/intent/follow?screen_name=immunomind)
[![CRAN](http://www.r-pkg.org/badges/version/immunarch?style=flat-square)](https://cran.r-project.org/package=immunarch)
[![Downloads_all](http://cranlogs.r-pkg.org/badges/grand-total/immunarch)](http://www.r-pkg.org/pkg/immunarch)
[![Downloads_week](http://cranlogs.r-pkg.org/badges/last-week/immunarch)](http://www.r-pkg.org/pkg/immunarch)
[![Issues](https://img.shields.io/github/issues/immunomind/immunarch?style=flat-square)](http://github.com/immunomind/immunarch/issues)
[![CI](https://gitlab.com/immunomind/immunarch/badges/master/pipeline.svg?style=flat-square)](https://gitlab.com/immunomind/immunarch/-/jobs)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3367200.svg)](https://doi.org/10.5281/zenodo.3367200)


# `immunarch` --- An R Package for Painless Bioinformatics Analysis of T-cell and B-cell Immune Repertoire Data


## Why `immunarch`?
- **Work with any type of data:** single-cell, bulk, data tables, databases --- you name it.
- **Community at the heart:** thrive in the community of almost 30,000 researchers and medical scientists worldwide, including researchers from **Pfizer, Novartis, Regeneron, Stanford, UCSF** and **MIT**.
- **One plot --- one line:** write a [whole PhD thesis in 8 lines of code](https://twitter.com/Nusob88/status/1127601201112129536) or reproduce almost any publication in 5-10 lines of `immunarch` code.
- **Be on the bleeding edge of science:** we regularly update `immunarch` with the latest methods. [Let us know what you need!](#help-the-community)
- **Automatic format detection and parsing** for all popular immunosequencing formats: from **MiXCR** and **ImmunoSEQ** to **10XGenomics** and **ArcherDX**.

And the **most important**: `immunarch` is not just a tool --- it is an **ecosystem**.


## From Berkeley with devotion

`immunarch` is brought to you by [ImmunoMind](https://immunomind.io) --- a [UC Berkeley SkyDeck](https://www.forbes.com/sites/avivalegatt/2019/01/07/launch-your-startup-at-these-five-college-incubators/) startup. 

We have been helping researchers to extract insights from sequencing data of T-cell and antibody repertoires since the inception of the RepSeq domain. Our bioinformatics tools are trusted by top universities including Stanford, UCSF, MIT, King's College London and big pharma companies including Pfizer and Novartis.


## Stay connected!

<form action="https://formspree.io/mjveopnl" method="POST">
  <br>
  <div class="group">      
      <input type="text" required name="_replyto" style="width:50%">
      <span class="highlight" style="width:50%; height:30%"></span>
      <span class="bar" style="width:50%"></span>
      <label>Leave your email and get the latest immunarch news</label>
      <br>
      <button class="bttn-fill bttn-md bttn-danger" type="submit">Connect!</button>
  </div>
  
  <input type="text" name="_gotcha" style="display:none" />
  <input type="hidden" name="_next" value="https://immunarch.com"/>
</form>

Feel free to follow us on [Twitter](https://twitter.com/immunomind) as well.

---

## Table of Contents

- [Introduction](#introduction)
- [Contact](#contact)
- [Installation](#installation)
- [Features](#features)
- [Quick Start](#quick-start)
- [Bugs and Issues](#bugs-and-issues)
- [Contribution](#help-the-community)
- [Citation](#citation)

## Introduction

`immunarch` is an R package designed to analyse T-cell receptor (TCR) and B-cell receptor (BCR) repertoires, aimed at medical scientists and bioinformaticians. The mission of `immunarch` is to make immune sequencing data analysis as effortless as possible and help you focus on research instead of coding. Follow us on [Twitter](https://twitter.com/immunomind) for news and updates.


## Contact
Create a ticket with a bug or question on [GitHub Issues](https://github.com/immunomind/immunarch/issues) to help the community help you and enrich it with your experience. If you need to send us a sensitive data, feel free to contact us via [support@immunomind.io](mailto:support@immunomind.io).


## Installation

### Latest release on CRAN
In order to install `immunarch` execute the following command:

```r
install.packages("immunarch")
```

That's it, you can start using `immunarch` now! See the [Quick Start](#quick-start) section below to dive into immune repertoire data analysis. If you run in any trouble with installation, take a look at the [Installation Troubleshooting](https://immunarch.com/articles/v1_introduction.html#installation-troubleshooting) section.

Note: there are quite a lot of dependencies to install with the package because it installs all the widely-used packages for data analysis and visualisation. You got both the AIRR data analysis framework and the full Data Science package ecosystem with only one command, making `immunarch` the entry-point for single-cell & immune repertoire Data Science.


### Latest release on GitHub
If the above command doesn't work for any reason, try installing `immunarch` directly from its repository:

```r
install.packages("devtools") # skip this if you already installed devtools
devtools::install_github("immunomind/immunarch")
```


### Latest pre-release on GitHub
Since releasing on CRAN is limited to one release per one-two months, you can install the latest pre-release version with bleeding edge features and optimisations directly from the code repository. In order to install the latest pre-release version, you need to execute only two commands:

```r
install.packages("devtools") # skip this if you already installed devtools
devtools::install_github("immunomind/immunarch", ref="develop")
```

You can find the list of releases of `immunarch` here: https://github.com/immunomind/immunarch/releases


## Features

1. Fast and easy manipulation of immune repertoire data:

    + The package automatically detects the format of your files---no more guessing what format is *that* file, just pass them to the package;
  
    + Supports all popular TCR and BCR analysis and post-analysis formats, including single-cell data: [ImmunoSEQ](https://www.adaptivebiotech.com/products-services/immunoseq/), [IMGT](http://www.imgt.org/IMGTindex/IMGTHighV-QUEST.php), [MiTCR](https://github.com/milaboratory/mitcr/), [MiXCR](https://milaboratory.com/software/mixcr/), [MiGEC](https://milaboratory.com/software/migec/), [MigMap](https://github.com/mikessh/migmap), [VDJtools](https://milaboratory.com/software/vdjtools/), [tcR](https://github.com/imminfo/tcr), [AIRR](http://docs.airr-community.org/en/latest/), [10XGenomics](https://support.10xgenomics.com/single-cell-vdj/datasets/), [ArcherDX](https://archerdx.com/immunology/). More coming in the future;

    + Works on any data source you are comfortable with: R data frames, data tables from [data.table](https://rdatatable.gitlab.io/data.table/), databases like [MonetDB](https://github.com/MonetDB), Apache Spark data frames via [sparklyr](https://spark.rstudio.com/);
    
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

### Use `immunarch` data

**1) Load the package and the data**

```r
library(immunarch)  # Load the package into R
data(immdata)  # Load the test dataset
```

**2) Calculate and visualise basic statistics**

```r
repExplore(immdata$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(immdata$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
```

**3) Explore and compare T-cell and B-cell repertoires**
```r
repOverlap(immdata$data) %>% vis()  # Build the heatmap of public clonotypes shared between repertoires
geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)  # Visualise the Chao1 diversity of repertoires, grouped by the patient status
```

### Use your own data

```r
library(immunarch)  # Load the package into R
immdata <- repLoad("path/to/your/data")  # Replace it with the path to your data. Immunarch automatically detects the file format.
```

### Advanced methods

For advanced methods such as clonotype annotation, clonotype tracking, kmer analysis and public repertoire analysis see "Tutorials".


## Bugs and Issues

The mission of `immunarch` is to make bulk and single-cell immune repertoires analysis painless. All bug reports, documentation improvements, enhancements and ideas are appreciated. Just let us know via [GitHub](https://github.com/immunomind/immunarch/issues) (preferably) or [support@immunomind.io](mailto:support@immunomind.io) (in case of private data).

Bug reports must: 

1. Include a short, self-contained R snippet reproducing the problem. 
2. Add a minimal data sample for us to reproduce the problem. In case of sensitive data you can send it to [support@immunomind.io](mailto:support@immunomind.io) instead of GitHub issues.
3. Explain why the current behavior is wrong/not desired and what you expect instead.
4. If the issue is about visualisations, please attach a picture to the issue. In other case we wouldn't be able to reproduce the bug and fix it.


## Help the community

Have an aspiration to help the community build the ecosystem of scRNAseq & AIRR analysis tools? Found a bug? A typo? Would like to improve a documentation, add a method or optimise an algorithm?

We are always open to contributions. There are two ways to contribute:

1. Create an issue [here](http://github.com/immunomind/immunarch/issues) and describe what would you like to improve or discuss.

2. Create an issue or find one [here](http://github.com/immunomind/immunarch/issues), fork the repository and make a pull request with the bugfix or improvement.


## Citation

ImmunoMind Team. (2019). immunarch: An R Package for Painless Bioinformatics Analysis of T-Cell and B-Cell Immune Repertoires. Zenodo. http://doi.org/10.5281/zenodo.3367200

BibTex:
```
@misc{immunomind_team_2019_3367200,
  author       = {{ImmunoMind Team}},
  title        = {{immunarch: An R Package for Painless Bioinformatics Analysis 
                    of T-Cell and B-Cell Immune Repertoires}},
  month        = aug,
  year         = 2019,
  doi          = {10.5281/zenodo.3367200},
  url          = {https://doi.org/10.5281/zenodo.3367200}
}
```

For EndNote citation import the [`immunarch-citation.xml`](https://gitlab.com/immunomind/immunarch/raw/master/immunarch-citation.xml?inline=false) file.

Preprint on BioArxiv is coming soon.


## License

The package is freely distributed under the AGPL v3 license. You can read more about it [here](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0)).

For commercial or server use, please contact ImmunoMind via [support@immunomind.io](mailto:support@immunomind.io) about solutions for biomarker data science of single-cell immune repertoires.

<!--
- Package modifications and feature implementations are issued promptly;
- Use *immunarch* team expertise in your projects;
- Priority email and call support;
- 100+ hours of consultations on the TCR & BCR repertoire analysis;
- Setup a cloud or cluster installation of *immunarch*, including the development of cloud *immunarch*-based software;
- If you need license other than the current, contact us.

Contact us at support@immunomind.io for more information.
-->

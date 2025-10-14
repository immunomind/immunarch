<div align="center">
  <h2>🧬 <code>immunarch</code> – <strong>Multi-modal immune repertoire analysis in R</strong></h2>
</div>

---

<div align="center">
  <a href="https://github.com/immunomind">
    <img alt="Ecosystem: ImmunoMind"
         src="https://img.shields.io/badge/ecosystem-ImmunoMind-orange?style=flat-square">
  </a>
  <a href="https://cran.r-project.org/package=immunarch">
    <img alt="CRAN Version"
         src="https://www.r-pkg.org/badges/version-ago/immunarch?style=flat-square">
  </a>
  <a href="https://www.r-pkg.org/pkg/immunarch">
    <img alt="CRAN Downloads (all time)"
         src="https://cranlogs.r-pkg.org/badges/grand-total/immunarch">
  </a>
  <a href="https://www.r-pkg.org/pkg/immunarch">
    <img alt="CRAN Downloads (last week)"
         src="https://cranlogs.r-pkg.org/badges/last-week/immunarch">
  </a>
  <a href="https://anaconda.org/conda-forge/r-immunarch">
    <img alt="Conda Version"
         src="https://anaconda.org/conda-forge/r-immunarch/badges/version.svg">
  </a>
  <a href="https://anaconda.org/conda-forge/r-immunarch">
    <img alt="Conda Total Downloads"
         src="https://anaconda.org/conda-forge/r-immunarch/badges/downloads.svg">
  </a>
  <a href="https://github.com/immunomind/immunarch/issues">
    <img alt="GitHub Issues"
         src="https://img.shields.io/github/issues/immunomind/immunarch?style=flat-square">
  </a>
  <a href="https://doi.org/10.5281/zenodo.3367200">
    <img alt="DOI"
         src="https://zenodo.org/badge/DOI/10.5281/zenodo.3367200.svg">
  </a>
</div>

<p align="center">
  <a href="https://immunomind.github.io/docs/tutorials/single_cell/">Tutorials</a>
  |
  <a href="https://immunomind.github.io/docs/api/reference/">API reference</a>
  |
  <a href=https://immunomind.github.io/docs/>Ecosystem</a>
  |
  Publication (coming soon...)
</p>

---

`immunarch` brings a comprehensive analytics toolkit to build reproducible analysis pipelines for Adaptive Immune Receptor Repertoire (AIRR) data with a particular focus on designing personalized immunotherapies and vaccines. Key features are:

- **Multi-modal immune profiling:** compute receptor- and repertoire-level statistics leveraging single-cell, spatial, immunogenicity or any other receptor annotations;

- **Immunomics at scale:** work seamlessly with datasets that don't fit in memory;

- **Immune biomarker discovery:** stratify cohorts and timepoints, derive repertoire signatures (diversity/clonality, V/J usage, similarity), and track antigen-annotated clonotypes;

- **Feature engineering:** build Machine Learning-ready feature tables (receptor-, ssample- and cohort-level) from core repertoire metrics and annotations, with consistent IDs/metadata for downstream statistics or modeling;

- **Modular, extendable, adaptable:** add new analyses and metrics via a extension API, and use adapters to interoperate with other AIRR tools and formats.


## 🤔 Why `immunarch`?

As immune repertoire sequencing becomes a mainstream technology, adopted by major platforms and integrated into more translational and clinical workflows, tooling expectations are changing rapidly.
The pace of innovation and data growth sometimes outstrips what even the most dedicated tool developers can deliver.

That's why it's the perfect moment to step back and rethink **how** and **why** we analyze AIRR data. 
Instead of racing to patch each new problem, we need to prepare for the **next epoch** of immunomics.

**What defines this next epoch?**

A massive shift in focus: from pure research towards biomarker discovery, personalized immunotherapies, and integration of immune repertoire data into real clinical decision-making.

Today's AIRR analysis must handle:

- **Multi-modal data:** bulk and single-cell V(D)J, spatial transcriptomics, gene expression, clinical metadata, and antigen specificity -- all together;

- **Massive scale:** experiments that move from gigabytes to tens or hundreds of gigabytes, or even terabytes;

- **Reproducibility and collaboration:** workflows that need to be shared, versioned, and rerun months or years later, sometimes by new teams.

With this new landscape, the **how** and **why** of AIRR data analysis are evolving:

- The focus is moving from "can I parse my data?" to "can I robustly extract insights, find biomarkers, and build ML-ready features for discovery or diagnostics?"

- It's not enough for a toolkit to just work. It needs to scale, interoperate, and empower new kinds of science, including complex Deep Learning and foundation models.

By taking a step back and rethinking the core "how" and "why" of AIRR analysis, `immunarch` prepares you for the next epoch of immunomics -- so your science is ready, whatever comes next.

---

> [!WARNING]
> `immunarch` is evolving towards `1.0` version and undergoing huge changes.
> Please check the updates here: https://github.com/immunomind/immunarch/issues/432
> 
> To install the latest pre-1.0 version, use `pak::pkg_install("immunomind/immunarch@0.9.1")`
>
> To install the latest 1.0 pre-release version, use `pak::pkg_install("immunomind/immunarch")`

---

- 🤔 [Why `immundata`?](#-why--immundata-)
- 📦 [Installation](#-installation)
- ⚡ [Quick Start](#-quick-start)
- [📄 Documentation](#-documentation)
- [🪲 Bugs and Issues](#-bugs-and-issues)
- 🏷 [About](#-about)
  - [Citation](#citation)
  - [License](#license)
  - [Author and contributors](#author-and-contributors)
  - [Commercial usage](#commercial-usage)

---

## 📦 Installation

### Prerequisites

Before installing any release or pre-release version of `immunarch`, please install `pak` that will simplify the installation of any package, not just `immunarch`:

```r
install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))
```

More info if needed is available on [pak website](https://pak.r-lib.org/#arrow_down-installation).

### Install the latest version

To install the latest release of `immunarch`, simply run:

```r
pak::pkg_install("immunomind/immunarch")
```

Mind that this will install the package from our GitHub instead of CRAN. This method is much preferred due to limitations of CRAN and reliance on other packages, which are distributed via `pak` as well.

### Other installation options

We will periodically release `immunarch` on CRAN. To install it from CRAN, run 

```r
pak::pkg_install("immunarch")
```

If you are willing to try unstable yet bleeding edge features, or if there are some hot fix for your open GitHub ticket, please install the development version:

```r
pak::pkg_install("immunomind/immunarch@dev")
```


## ⚡ Quick Start

```r
# Install `pak` - a blazingly-fast package manager
install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))

# Install and load immunarch along with the pre-packaged data
pak::pkg_install("immundata/immunarch")
library(immunarch)
immdata <- get_test_idata() |> agg_repertoires("Therapy")

# In just 4 lines of code you get the essential AIRR statistics
airr_stats_genes(immdata, gene_col = "v_call") |> vis()
airr_public_jaccard(immdata) |> vis()
airr_diversity_pielou(immdata) |> vis()
airr_diversity_chao1(immdata) |> vis()
airr_clonality_prop(immdata)

# Use your own data by reading sample files from the metadata file
mdtable <- read_metadata("data/metadata.csv")
immdata <- read_repertoires("<metadata>", metadata = mdtable)

# Use your own data by reading sample files directly
mdtable <- read_metadata("data/metadata.csv")
immdata <- read_repertoires("data/*.tsv.gz", metadata = mdtable)
```

Oh, and one small thing. Even if you have tens of gigabytes of the data, you won't need to adapt the code to a server.
The code will be exactly the same — `immunarch` got your back thanks to [`immundata`](https://github.com/immunomind/immundata/).


## 📄 Documentation

To get a list of available methods and their descriptions, run the default help command in R on specific functions or on a function family prefix:

```r
# This is the same
?airr_stats
# as this
?airr_stats_genes

# Basic statistics - gene usage, length distribution
?airr_stats

# Public receptor indices - overlap, jaccard, morisita
?airr_public

# Clonality analysis - clonal lines, occupied space, clonal space homeostasis
?airr_clonality

# Diversity analysis - pielou, shannon, chao1
?airr_diversity

# ... more to come ...
```

More detailed documentation, guides and comprehensive tutorials are available on the ecosystem website: [https://immunomind.github.io/docs/](https://immunomind.github.io/docs/).


## 🪲 Bugs and Issues

The mission of `immunarch` is to make bulk and single-cell immune repertoires analysis painless. All bug reports, documentation improvements, enhancements and ideas are appreciated. Just let us know via [GitHub](https://github.com/immunomind/immunarch/issues) (preferably) or [support@immunomind.com](mailto:support@immunomind.com) (in case of private data).

Bug reports must: 

1. Include a short, self-contained R snippet reproducing the problem. 
2. Add a minimal data sample for us to reproduce the problem. In case of sensitive data you can send it to [support@immunomind.com](mailto:support@immunomind.com) instead of GitHub issues.
3. Explain why the current behavior is wrong/not desired and what you expect instead.
4. If the issue is about visualisations, please attach a picture to the issue. In other case we wouldn't be able to reproduce the bug and fix it.

We are always open to contributions. There are three ways to contribute:

1. Create an issue [here](https://github.com/immunomind/immunarch/issues) and describe what would you like to improve or discuss.

2. Create an issue or find one [here](https://github.com/immunomind/immunarch/issues), fork the repository and make a pull request with the bugfix or improvement.

3. Find an existing issue and help others resolve this.


## 🏷 About

### Citation

> Temporary citation is below. The main manuscript is in preparation. Preprint on BioArxiv is coming soon as of 2025.

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


### License

The package is freely distributed under the Apache-2.0 license. You can read more about it [here](https://www.tldrlegal.com/license/apache-license-2-0-apache-2-0).

### Author and contributors 

- **Vadim I. Nazarov – main author and developer**
- Vasily Tsvetkov
- Aleksandr Popov
- Ivan Balashov

### Commercial usage 

`immunarch` is free to use for commercial usage as per Apache-2.0 license. However, corporate users will not get a prioritized support for `immunarch`- or AIRR-related issues. The priority of open-source tool `immunarch` is open-source science.

If you are looking for prioritized support and setting up your data pipelines, consider contacting [Vadim Nazarov](https://www.linkedin.com/in/vdnaz/) for commercial consulting / support options / workshops and training sessions / designing data platforms and machine learning systems for multi-omics / or anything related.

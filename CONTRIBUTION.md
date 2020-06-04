# Contributing to immunarch

This document outlines how to propose a change to immunarch.

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitLab web IDE functionality / GitHub web interface, so long as the changes are made in the _source_ file. Please make sure to create a Pull Request instead of commiting directly to the `master` branch.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it's a problem. If you've found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex). If you need to send a data and can't do
that via GitHub, please note it in the issue.

### Version naming

We employ Major.Minor.Patch. E.g., 0.4.0.

Hotfixes: Major.Minor.Patch → Major.Minor.(Patch + 1)

Guidelines: [http://r-pkgs.had.co.nz/release.html](http://r-pkgs.had.co.nz/release.html)

### Commit naming

We follow the [following guidelines for commit naming](https://medium.com/@kevinkreuzer/the-way-to-fully-automated-releases-in-open-source-projects-44c015f38fd6).

#### Commit types

There are eight types of commits: chore, docs, feat, fix, refactor, test, perf, style.

#### Commit scopes

Chanes in analysis and visualisation functions: diversity, overlap, public-repertoire, clonality, gene-usage, explore, kmers, spectratype, dynamics, tools

Changes in parsing: io

Changes in databases support: db

Changes in additional functions such as general statistics functions: utility

Changes in NAMESPACE, DESCRIPTION, citations, ISSUE_TEMPLATE.md, etc., without README: upkeep

Changes in README and vignettes: vignette

Changes in Continuous Integration: ci

Changes in Shiny applications: shiny

Changes in the website: web

### How to create an Issue

We have a rich list of templates for Issues [here](https://github.com/immunomind/immunarch/tree/master/.github/ISSUE_TEMPLATE).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the CI (Continuous Integration) build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow our style guide that is the tidyverse [style guide](http://style.tidyverse.org) You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2) for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the current
development version header describing the changes made followed by your GitHub
username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](code-of-conduct.md). By participating in this project you agree to
abide by its terms.

(adapted from [dplyr](https://github.com/tidyverse/dplyr/blob/master/.github/CONTRIBUTING.md) contribution guide)

# Contributing to immunarch

This outlines how to propose a change to immunarch.

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitLab web IDE functionality / GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it's a problem. If you've found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex). If you need to send a data and can't do
that via GitHub, please note it in the issue.

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the CI (Continuous Integration) build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow our style guide that is the tidyverse [style guide](http://style.tidyverse.org)
with a very minor addition of using `=` instead of `<-` in variable assignment.
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
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

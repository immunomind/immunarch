# immunarch 0.4.3

## Documentation

Make README more clear with the focus on installation.

Add a vignette for clonotype tracking.

Add a vignette for kmer analysis.

## Work in progress

Kmers - analysis + plots, seqlogo / textlogo

`pubRepStatistics` function for computation of public clonotypes across groups and repertoires.

`vis_public_frequencies`, `vis_public_clonotypes`, `vis_upset` for public repertoire visualisations.


# immunarch 0.4.1 - 0.4.2

Minor releases with bug fixes and minor improvements.


# immunarch 0.4.0

## Features

The `trackClonotype` function for tracking clonotypes.

Visualisations for `trackClonotype`.

Speed up Morisita-Horn index computation for up to 3 times.

Add the incremental overlap function `inc_overlap`, available from `repOverlap`.

Add the downsampling option to incremental overlaps.

Add visualisations for incvremental overlaps.

Remove the `prop_sample` function completely.

Move all downsampling, resampling and sampling procedures to the `repSample` function.

Add normalisation to rarafection.

Remove the .dup argument from `repOverlap`, all equal clonotypes are now always merged and their counts summed up.

Remove .quant from `repOverlap`, the column is now automatically detected.

Add an argument .add.layer to `vis_hist` to add any additional ggplot2 layers to the plots in the output grid.

Add grouping in `vis_hist` if data is grouped and `.grid` is T.

Update MiXCR parser (yet again) to read MiXCR files from the May 2019 release.

Add `.target` argument to visualisation of incremental overlaps.

Add `.clones` to `repExplore`.

Add more MiXCR file variants parsing.

`repLoad` returns sorted by "Clones" data frames now.

## Bug fixes

Fix `geneUsage` when the first two columns with gene usages were swapped.

Fix `repExplore` doesn't work with a single repertoire.

Fix `vis_heatmap` doesn't work with `geneUsage` output.

Fix a bug in computation of Morisita-Horn index.

Fix `vis_bar` and `vis_box` doesn't work with numeric grouping variables.

Fix `repOverlap` failing when working on data tables with `morisita`.

Fix `repOverlap` failing to work with data tables on `public`, `overlap`, `jaccard` and `tversky`.

Fix `repLoad` failing when parsing MiXCR files with zero clonotypes.

Fix a bug in incorrect grouping in visualisations.

Fix a bug when clonal homeostasis and clonotype tracking don't work properly with filtered coding.

Fix D50, "top" and "clonal.prop" from repClonality returns wrong values when the input data frame is not sorted.

## Minor updates

Remove the "wei" option from `geneUsage` becase it's useless.

New dependency `ggseqlogo` for visualising of seq-logo plots.

Add `.transpose` to `vis_heatmap`.

`repOverlap` default .col value is "aa" for comfortable usage.

Remove "fill" aesthetics warnings from `vis_heatmap`.

Remove warnings "In parse_fun(.path[i]) : NAs introduced by coercion" when parsing MiXCR files.

A lot or minor fixing and documentation improving to prepare for the CRAN release.

Remove short function names.


# immunarch 0.3.3

## Features

Parser for ArcherDX.

Update a parser for MiXCR to make it work with the "targetSequences" column format.

Update the coding function family to make it work with CDR3 amino acid sequences only.

## Minor updates

Replace all "unresolved" genes in ImmunoSEQ parsed files with NAs.

Add an argument for the color palette to `vis_heatmap2`.

Make `.grid=F` by default in `vis_hist` for gene usage analysis.

Fix parsing functions to not remove strings after dots in filenames.

Update the post-parsing processing subroutine to remove all characters except for amino acid alphabet, `*` and `~` for compatability with all `immunarch`'s functions.

## Documentation

Add a documentation to the coding function family.

Update README.


# immunarch 0.3.2

Remove MonetDBLite from dependencies because it got removed from CRAN.

## Bug fixes
* Fix a bug in MiXCR parser.


# immunarch 0.3.1

## Features
* Boxplots and barplots are now support statistical tests via the `.test` argument.

* Add parsers for old VDJtools formats.

## Documentation
* Update docs and vignettes with statistical tests information.

* Add a note for list names to vignettes.

* Documentation for clustering.

* Documentation for dimension reduction.

* Minor fix for `repOverlap` documentation.

## Bug fixes
* Fix a grouping bug in visualisations.

* Fix statistical tests from `ggpubr`.

* Fix for `geneUsage` with `.type="family"`.


# immunarch 0.3.0

## Features
* `fixVis` now supports the following legends: size, shape, color, fill, linetype.

* `fixVis` can plot figures to R console / RStudio "Plots" tab.

* `fixVis` now supports the number of columns in legends.

* Support for the AIRR file format.

* Experimental support for the 10xGenomics format.

* Save and load `immunarch` format via `repSave` and `repLoad`.

* Save and load VDJtools format via `repSave` and `repLoad`.

## Bug fixes
* `.a` and `.b` didn't passed to Tversky index.

* `fixVis` - fix a bug when users apply X/Y settings to the other axis.

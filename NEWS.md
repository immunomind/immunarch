# Releases

---

## immunarch 0.11.0 (current development version)

This release adds receptor tracking, publicness analysis, sequence clustering,
and immunoglobulin mutation analysis to the `ImmunData` API. It also expands
repertoire similarity and adds rarefaction and extrapolation.

### Breaking changes

* Renamed `airr_stats_chains()`, `airr_stats_genes()`, and `airr_stats_lengths()`
  to `airr_desc_chains()`, `airr_desc_genes()`, and `airr_desc_lengths()`.
* Renamed `airr_public_intersection()` and `airr_public_jaccard()` to
  `repsim_intersection()` and `repsim_jaccard()`.
* If you used `clust_cc()` in an earlier development build, pass the contracted
  graph from `dist_hamm(..., mode = "graph")` as `edges`.

### New features and improvements

* Added `vdj_track()` to follow selected receptors across samples, time points,
  or tissues. Use `output = "df"` to get a table with zero counts for receptors
  absent from a repertoire. Plot the tracked abundances with `vis()`.
* Added `annotate_public()` to measure receptor incidence and abundance across
  repertoires and strata. Use `vdj_public()` and `make_public_filter()` to select
  receptors by incidence, count, or proportion. Use `vis()` to plot incidence
  distributions and public receptor selections.
* Added `dist_hamm()` for grouped Hamming distances and `clust_cc()` for
  connected-component clustering. Distance calculations support thresholds,
  sampling within groups, and separate comparisons by locus for paired chains.
* Added Chao-Jaccard, Morisita-Horn, and Bray-Curtis repertoire similarity.
  Shared-receptor counts and Jaccard similarity skip unshared receptors
  when building repertoire pairs.
  All five similarity methods support pairwise dot plots with `vis()`.
* Updated descriptive statistics for paired chains. Added `vis()` plots for
  chain counts and for rank-abundance, rank bins, and clonal space.

### Receptor annotation

* Added `annotate_with_database()` to attach a reference database label by
  matching selected receptor features.
* Added `annotate_ig_clones()` to assign immunoglobulin clones from heavy-chain
  or paired-chain data. Use `read_germline_references()` to load IMGT references
  and `annotate_ig_germlines()` to reconstruct clone-level germlines.
* Added `annotate_ig_shm()` to count replacement and silent mutations against
  aligned germlines. Use `airr_ig_shm()` to summarize mutation counts and
  frequencies by repertoire, with summaries across sequences and across clones.

### Rarefaction and extrapolation

* Added `airr_diversity_rarefaction()` and support for plotting its results with
  `vis()`. The function requires positive integer clonotype counts.
* Interpolation uses the exact hypergeometric expectation for sampling without
  replacement. This corrects underestimated richness in earlier development
  builds.
* Extrapolation uses the finite-sample Chao estimator and defaults to twice
  the observed sample size. Confidence bounds use an abundance-distribution
  bootstrap and both probabilities supplied through `quantile`.

### Bug fixes

* Fixed Chao1 richness estimates for repertoires with doubletons and no
  singletons. `airr_diversity_dxx(..., perc = 100)` now returns observed
  richness even when cumulative proportions round below 100%.
* Repertoire similarity counts each receptor once per repertoire when
  calculating shared receptor counts or Chao-Jaccard abundances.
* AIRR methods now honor `options(immunarch.autojoin = FALSE)` and avoid
  duplicating repertoire metadata columns already present in the output.
* Fixed DuckDB 2.0 compatibility for Hamming distances, clone annotation,
  and somatic hypermutation annotation.

### Deprecations

* The January 2027 release, immunarch 0.12.0, will remove legacy functions
  that currently raise deprecation warnings. Run
  `immunarch::get_immunarch_news()` for migration help.
  [Migration guide](https://immunomind.github.io/docs/tutorials/migration/),
  [single-cell tutorial](https://immunomind.github.io/docs/tutorials/single_cell/)
* Include the warning and a small reproducible example when reporting a
  migration problem. [GitHub issues](https://github.com/immunomind/immunarch/issues)

### Documentation and maintenance

* Raised the minimum versions to immundata 0.1.0, dplyr 1.2.0, and duckplyr
  1.2.0. Removed Rcpp and reshape2, and moved more legacy dependencies to
  `Suggests`.
* Added examples and tests for paired-chain summaries, repertoire similarity,
  receptor tracking, clustering, and immunoglobulin annotation.

---

## immunarch 0.10.3 (October 2025)

This release previews the 1.0 API. It introduces `ImmunData` workflows and
keeps legacy functions available during migration.

### Breaking changes

* Raised the minimum R version from 4.0.0 to 4.2.0.
* Removed the bundled `scdata` example. The single-cell tutorial provides data
  and examples for the `ImmunData` workflow.
  [Single-cell tutorial](https://immunomind.github.io/docs/tutorials/single_cell/)

### New features and improvements

* Added analyses for `ImmunData` from immundata. The duckplyr and DuckDB
  backend supports datasets larger than memory. Use
  `immundata::read_repertoires()` to load bulk or single-cell data, or
  `immundata::from_immunarch()` to convert results from `repLoad()`.
* Added `airr_stats_*()` for chain counts, sequence lengths, and gene usage.
  The new data model supports paired chains and links receptors to single-cell
  or spatial metadata.
* Added `airr_clonality_*()` for rank-abundance, rank bins, and clonal space.
  Added `annotate_clonality_rank()` and `annotate_clonality_prop()` to label
  receptors by rank or abundance proportion.
* Added `airr_diversity_*()` for richness, entropy, evenness, and Hill numbers.
  Added `airr_public_intersection()` and `airr_public_jaccard()` for pairwise
  repertoire overlap.
* Extended `vis()` to selected results from the new API, including diversity
  measures and pairwise repertoire overlap.

### Documentation and maintenance

* Marked legacy functions such as `repLoad()` and `repDiversity()` as deprecated
  in their help pages. Added `get_immunarch_news()` and
  `list_immunarch_news()` to share migration updates.
  [Migration guide](https://immunomind.github.io/docs/tutorials/migration/)
* Added vignettes on the new API and DuckDB thread settings. Moved several
  dependencies for legacy plotting, clustering, and dimensionality reduction
  to `Suggests`.

---

## immunarch 0.9.1 (March 2024)

This is the last release based solely on the legacy API, including `repLoad()`,
`repDiversity()`, and `repOverlap()`. Install it with
`pak::pak("immunomind/immunarch@0.9.1")`.

### Bug fixes

* Fixed MiXCR import for files without a `cloneId` column. The parser assigns
  identifiers from row numbers and detects non-tab delimiters.
* Fixed AIRR import for files missing optional CDR and framework region
  columns. The parser fills the missing columns with `NA`.
* Fixed sample name extraction from paths with Windows separators.
  `repLoad()` now separates mismatched sample names with commas when
  reporting metadata problems.
* Registered S3 methods for the overlap coefficient, Jaccard and Tversky
  indices, and cosine similarity. Updated `vis_box()` and `vis_bar()` to use
  ggplot2's `after_stat()` for adjusted p-value labels.

### Documentation and maintenance

* Updated installation and data-loading examples. Fixed broken links in the
  README and vignettes.

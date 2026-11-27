

# Build and annotate immunoglobulin clonal families

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/#L)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

Build clonal families from B-cell receptors and assign them to
receptors.

Receptors are compared within a subject or repertoire. A heavy-chain
edge requires a shared V/J union group, equal junction length, and
normalized junction Hamming distance at or below
<code>heavy_threshold</code>. Paired mode independently builds the
corresponding light-chain graph and keeps receptors together only when
they share both heavy- and light-chain families.

## Usage

<pre><code class='language-R'>annotate_ig_clones(
  idata,
  heavy_threshold,
  light_threshold = NULL,
  mode = c("auto", "heavy", "paired"),
  by_scope = "&lt;auto&gt;",
  junction_col = "junction",
  v_col = "v_call",
  j_col = "j_call",
  heavy_locus = "IGH",
  clone_col = "clone_id"
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>immundata::ImmunData. Paired or heavy-chain IG receptor data.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="heavy_threshold">heavy_threshold</code>, <code
id="light_threshold">light_threshold</code></td>
<td>Numeric scalar or data frame. Normalized Hamming-distance threshold
strictly between <code>0</code> and <code>1</code>, or a complete table
containing the <code>by_scope</code> column and a numeric
<code>threshold</code> column. <code>light_threshold</code> defaults to
<code>heavy_threshold</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="mode">mode</code></td>
<td>Character scalar. Annotation strategy: <code>"heavy"</code> uses
only the heavy chain, <code>"paired"</code> combines independently
inferred heavy- and light-chain families, and <code>"auto"</code>
selects paired mode for a paired IG receptor schema and heavy mode
otherwise.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by_scope">by_scope</code></td>
<td>Character scalar. Column that prevents clonal families from crossing
biological subjects or repertoires. <code>"&lt;auto&gt;"</code> uses the
first available subject-like column and otherwise the repertoire
identifier.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="junction_col">junction_col</code>, <code id="v_col">v_col</code>,
<code id="j_col">j_col</code></td>
<td>Character scalar. Columns containing the nucleotide junction, V
call, and J call. V and J calls are grouped at gene level, allowing
overlapping ambiguous assignments.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="heavy_locus">heavy_locus</code></td>
<td>Character scalar. Value identifying the heavy-chain locus. The
default is <code>"IGH"</code>; all other loci declared by the receptor
schema are treated as light chains in paired mode.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="clone_col">clone_col</code></td>
<td>Character scalar. Name of the combined clone annotation. The default
is <code>"clone_id"</code>.</td>
</tr>
</tbody>
</table>

## Details

Connected components use single linkage: when A connects to B and B
connects to C, all three receptors belong to one chain family even if A
and C are farther apart than the threshold. Hamming distance compares
only junctions of equal length, so it does not represent insertions or
deletions.

Junction nucleotides are converted to uppercase. Junctions containing
characters other than A, C, G, and T cannot create edges and remain
isolated.

For grouping, V/J calls are split on commas or <code>|</code>,
whitespace and allele suffixes are removed, and duplicate calls are
ignored. Calls sharing both a V candidate and a J candidate are
connected transitively within each scope and locus, before splitting by
junction length or applying distance thresholds. Missing or empty gene
calls cannot connect groups. Original calls and receptor IDs are
preserved. The second example separately shows preprocessing when
normalized calls should also define receptor identity.

## Value

immundata::ImmunData. The original rows and columns plus
<code>clone_id</code> (or the name supplied through
<code>clone_col</code>), <code>heavy_family_id</code>, and
<code>light_family_id</code>. Both chain rows of a paired receptor
receive the same labels. In heavy mode, <code>clone_id</code> equals
<code>heavy_family_id</code> and <code>light_family_id</code> is
missing.

## Examples

``` r
library("immunarch")




# Limit the number of threads used by the underlying DB for this example.
# Generally, you should NOT do this in your session.
db_exec("SET threads TO 1")

# Six chain observations representing three paired B-cell receptors.
ig_chains <- tibble::tribble(
  ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
  "S1", "c1", "P1", "IGH", "AAAA", "IGHV1", "IGHJ1", 10L,
  "S1", "c1", "P1", "IGK", "CCCC", "IGKV1", "IGKJ1", 8L,
  "S1", "c2", "P1", "IGH", "AAAT", "IGHV1", "IGHJ1", 10L,
  "S1", "c2", "P1", "IGK", "CCCT", "IGKV1", "IGKJ1", 8L,
  "S1", "c3", "P1", "IGH", "GGGG", "IGHV2", "IGHJ2", 10L,
  "S1", "c3", "P1", "IGL", "TTTT", "IGLV1", "IGLJ1", 8L
)
input_file <- tempfile(fileext = ".tsv")
output_folder <- tempfile("ig-clones-")
readr::write_tsv(ig_chains, input_file)

ig_schema <- make_receptor_schema(
  features = c("junction", "v_call", "j_call"),
  chains = c("IGH", "IGK|IGL")
)

idata <- read_repertoires(
  path = input_file,
  schema = ig_schema,
  barcode_col = "cell",
  locus_col = "locus",
  umi_col = "umi",
  repertoire_schema = "sample",
  rename_columns = NULL,
  output_folder = output_folder,
  verbose = FALSE
)

# c1 and c2 form one family; c3 remains separate.
annotated <- annotate_ig_clones(idata, heavy_threshold = 0.25)

annotated
```

    #> # A duckplyr data frame: 9 variables
    #>   imd_receptor_id junction.x v_call.x j_call.x locus.x junction.y v_call.y
    #>             <int> <chr>      <chr>    <chr>    <chr>   <chr>      <chr>   
    #> 1               3 GGGG       IGHV2    IGHJ2    IGH     TTTT       IGLV1   
    #> 2               1 AAAA       IGHV1    IGHJ1    IGH     CCCC       IGKV1   
    #> 3               2 AAAT       IGHV1    IGHJ1    IGH     CCCT       IGKV1   
    #> # ℹ 2 more variables: j_call.y <chr>, locus.y <chr>

    #> # A duckplyr data frame: 20 variables
    #>   imd_receptor_id imd_chain_id sample cell  subject locus junction v_call j_call
    #>             <int>        <int> <chr>  <chr> <chr>   <chr> <chr>    <chr>  <chr> 
    #> 1               1            1 S1     c1    P1      IGH   AAAA     IGHV1  IGHJ1 
    #> 2               1            2 S1     c1    P1      IGK   CCCC     IGKV1  IGKJ1 
    #> 3               2            3 S1     c2    P1      IGH   AAAT     IGHV1  IGHJ1 
    #> 4               2            4 S1     c2    P1      IGK   CCCT     IGKV1  IGKJ1 
    #> 5               3            5 S1     c3    P1      IGH   GGGG     IGHV2  IGHJ2 
    #> 6               3            6 S1     c3    P1      IGL   TTTT     IGLV1  IGLJ1 
    #> # ℹ 11 more variables: umi <dbl>, imd_filename <chr>, imd_barcode <chr>,
    #> #   imd_n_chains <dbl>, imd_count <dbl>, imd_repertoire_id <int>,
    #> #   imd_proportion <dbl>, n_repertoires <int>, clone_id <int>,
    #> #   heavy_family_id <int>, light_family_id <int>

    #> # A tibble: 1 × 4
    #>   imd_repertoire_id sample n_barcodes n_receptors
    #>               <int> <chr>       <int>       <int>
    #> 1                 1 S1              3           3

``` r
# Use heavy-chain similarity alone when the scientific question requires it.
annotate_ig_clones(
  idata,
  heavy_threshold = 0.25,
  mode = "heavy"
)
```

    #> # A duckplyr data frame: 9 variables
    #>   imd_receptor_id junction.x v_call.x j_call.x locus.x junction.y v_call.y
    #>             <int> <chr>      <chr>    <chr>    <chr>   <chr>      <chr>   
    #> 1               3 GGGG       IGHV2    IGHJ2    IGH     TTTT       IGLV1   
    #> 2               1 AAAA       IGHV1    IGHJ1    IGH     CCCC       IGKV1   
    #> 3               2 AAAT       IGHV1    IGHJ1    IGH     CCCT       IGKV1   
    #> # ℹ 2 more variables: j_call.y <chr>, locus.y <chr>

    #> # A duckplyr data frame: 20 variables
    #>   imd_receptor_id imd_chain_id sample cell  subject locus junction v_call j_call
    #>             <int>        <int> <chr>  <chr> <chr>   <chr> <chr>    <chr>  <chr> 
    #> 1               1            1 S1     c1    P1      IGH   AAAA     IGHV1  IGHJ1 
    #> 2               1            2 S1     c1    P1      IGK   CCCC     IGKV1  IGKJ1 
    #> 3               2            3 S1     c2    P1      IGH   AAAT     IGHV1  IGHJ1 
    #> 4               2            4 S1     c2    P1      IGK   CCCT     IGKV1  IGKJ1 
    #> 5               3            5 S1     c3    P1      IGH   GGGG     IGHV2  IGHJ2 
    #> 6               3            6 S1     c3    P1      IGL   TTTT     IGLV1  IGLJ1 
    #> # ℹ 11 more variables: umi <dbl>, imd_filename <chr>, imd_barcode <chr>,
    #> #   imd_n_chains <dbl>, imd_count <dbl>, imd_repertoire_id <int>,
    #> #   imd_proportion <dbl>, n_repertoires <int>, clone_id <int>,
    #> #   heavy_family_id <int>, light_family_id <int>

    #> # A tibble: 1 × 4
    #>   imd_repertoire_id sample n_barcodes n_receptors
    #>               <int> <chr>       <int>       <int>
    #> 1                 1 S1              3           3

``` r
# Canonicalize ambiguous V/J calls lazily before receptors are constructed.
# For example, "IGHV2*01, IGHV1*02" becomes "IGHV1|IGHV2".
canonicalize_vj <- function(data, ...) {
  data |>
    dplyr::mutate(dplyr::across(
      c(v_call, j_call),
      ~ dd$array_to_string(
        dd$list_sort(dd$list_distinct(dd$string_split(
          dd$regexp_replace(
            dd$regexp_replace(.x, "\\*[0-9]+", "", "g"),
            "\\s+", "", "g"
          ),
          ","
        ))),
        "|"
      )
    ))
}

messy_chains <- ig_chains |>
  dplyr::mutate(
    v_call = dplyr::case_when(
      cell == "c2" & locus == "IGH" ~ " IGHV1*02, IGHV1*01 ",
      TRUE ~ paste0(v_call, "*01")
    ),
    j_call = paste0(j_call, "*01")
  )
messy_file <- tempfile(fileext = ".tsv")
canonical_folder <- tempfile("ig-clones-canonical-")
readr::write_tsv(messy_chains, messy_file)

canonical_idata <- read_repertoires(
  path = messy_file,
  schema = ig_schema,
  barcode_col = "cell",
  locus_col = "locus",
  umi_col = "umi",
  repertoire_schema = "sample",
  rename_columns = NULL,
  preprocess = list(canonicalize_vj = canonicalize_vj),
  output_folder = canonical_folder,
  verbose = FALSE
)

annotate_ig_clones(canonical_idata, heavy_threshold = 0.25)
```

    #> # A duckplyr data frame: 9 variables
    #>   imd_receptor_id junction.x v_call.x j_call.x locus.x junction.y v_call.y
    #>             <int> <chr>      <chr>    <chr>    <chr>   <chr>      <chr>   
    #> 1               3 GGGG       IGHV2    IGHJ2    IGH     TTTT       IGLV1   
    #> 2               1 AAAA       IGHV1    IGHJ1    IGH     CCCC       IGKV1   
    #> 3               2 AAAT       IGHV1    IGHJ1    IGH     CCCT       IGKV1   
    #> # ℹ 2 more variables: j_call.y <chr>, locus.y <chr>

    #> # A duckplyr data frame: 20 variables
    #>   imd_receptor_id imd_chain_id sample cell  subject locus junction v_call j_call
    #>             <int>        <int> <chr>  <chr> <chr>   <chr> <chr>    <chr>  <chr> 
    #> 1               1            1 S1     c1    P1      IGH   AAAA     IGHV1  IGHJ1 
    #> 2               1            2 S1     c1    P1      IGK   CCCC     IGKV1  IGKJ1 
    #> 3               2            3 S1     c2    P1      IGH   AAAT     IGHV1  IGHJ1 
    #> 4               2            4 S1     c2    P1      IGK   CCCT     IGKV1  IGKJ1 
    #> 5               3            5 S1     c3    P1      IGH   GGGG     IGHV2  IGHJ2 
    #> 6               3            6 S1     c3    P1      IGL   TTTT     IGLV1  IGLJ1 
    #> # ℹ 11 more variables: umi <dbl>, imd_filename <chr>, imd_barcode <chr>,
    #> #   imd_n_chains <dbl>, imd_count <dbl>, imd_repertoire_id <int>,
    #> #   imd_proportion <dbl>, n_repertoires <int>, clone_id <int>,
    #> #   heavy_family_id <int>, light_family_id <int>

    #> # A tibble: 1 × 4
    #>   imd_repertoire_id sample n_barcodes n_receptors
    #>               <int> <chr>       <int>       <int>
    #> 1                 1 S1              3           3

``` r
unlink(messy_file)
unlink(canonical_folder, recursive = TRUE)

unlink(input_file)
unlink(output_folder, recursive = TRUE)
```

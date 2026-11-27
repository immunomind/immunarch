

# Connected-component clustering for receptor graphs

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

A family of functions for clustering sparse receptor graphs represented
by an <code>ImmunData</code> object and an
<code>immunarch_sequence_graph</code>.

<code>clust_cc()</code> finds receptor-level connected components from
the contracted sequence graph returned by <code>dist_hamm()</code> with
<code>mode = "graph"</code>.

## Usage

<pre><code class='language-R'>clust_cc(
  idata,
  edges,
  by = NULL,
  cluster_col = "cluster_id",
  size_col = "cluster_size"
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>An immundata::ImmunData object. Its receptor annotations provide the
graph nodes through the canonical receptor-ID column and, when
requested, repertoire metadata.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="edges">edges</code></td>
<td>An <code>immunarch_sequence_graph</code> returned by
<code>dist_hamm()</code> with <code>mode = "graph"</code>. Ordinary
receptor edge tables are not accepted.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>A character vector of column names or <code>NULL</code>. These
columns come from the receptor annotations in <code>idata</code>; each
distinct combination is an independent clustering scope. The scope may
be broader than the comparison groups used by <code>dist_hamm()</code>.
The canonical receptor-ID column cannot also be used in
<code>by</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="cluster_col">cluster_col</code></td>
<td>A single character string giving the output component-identifier
column name. Within each <code>by</code> group, the identifier is the
minimum original receptor identifier in that component. Defaults to
<code>"cluster_id"</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="size_col">size_col</code></td>
<td>A single character string giving the output component-size column
name. The size is the number of distinct original receptors in the
component. Defaults to <code>"cluster_size"</code>.</td>
</tr>
</tbody>
</table>

## Details

Sequence-node distances are applied only where both nodes participate in
the same clustering scope. If one receptor belongs to several sequence
nodes in a scope, those nodes are joined by identity edges before
components are computed. Receptors without a usable sequence remain
singleton components.

## Value

A lazy duckplyr table with one row per distinct, non-missing receptor ID
per clustering scope and the following columns:

<ul>
<li> the grouping columns supplied through <code>by</code>, in the same order;

</li>
<li> the canonical receptor identifier;

</li>
<li> the component identifier, named by <code>cluster_col</code>; and

</li>
<li> the number of distinct receptors in that component, named by <code>size_col</code>.

</li></ul>

Call <code>dplyr::collect()</code> when an in-memory tibble is needed.

## Examples

``` r
library("immunarch")


# Limit the number of threads used by the underlying DB for this example.
# Generally, you should NOT do this in your session.
db_exec("SET threads TO 1")

idata <- immundata::get_test_idata()

# Connect receptors whose CDR3 amino-acid sequences differ by at most one
# residue. Request the contracted graph used by clust_cc().
graph <- dist_hamm(idata, max_dist = 1, mode = "graph")
res <- clust_cc(
  idata = idata,
  edges = graph
)
res
```

    #> # A duckplyr data frame: 3 variables
    #>    imd_receptor_id cluster_id cluster_size
    #>              <int>      <int>        <dbl>
    #>  1             883        883            1
    #>  2             625        625            1
    #>  3            1000       1000            1
    #>  4            1151       1151            1
    #>  5            1540       1540            1
    #>  6             812        812            1
    #>  7            1595       1595            1
    #>  8             301        301            1
    #>  9            1175       1175            1
    #> 10            1326       1326            1
    #> # ℹ more rows

``` r
res <- res |> collect()


# Build independent graphs within each repertoire and V-gene group. The
# same columns must be supplied to by in both functions.
group_cols <- c("imd_repertoire_id", "v_call")
graph_by_group <- dist_hamm(
  idata,
  by = group_cols,
  min_sim = 0.8,
  mode = "graph"
)

res_by_group <- clust_cc(
  idata = idata,
  edges = graph_by_group,
  by = group_cols
)

res_by_group |>
  dplyr::filter(cluster_size > 1L) |>
  dplyr::arrange(dplyr::desc(cluster_size)) |>
  dplyr::collect()
```

    #> # A tibble: 57 × 5
    #>    imd_repertoire_id v_call      imd_receptor_id cluster_id cluster_size
    #>  *             <int> <chr>                 <int>      <int>        <dbl>
    #>  1                 1 TRBV20-1*01             311         49            4
    #>  2                 1 TRBV20-1*01             313         49            4
    #>  3                 1 TRBV20-1*01             588         49            4
    #>  4                 1 TRBV20-1*01              49         49            4
    #>  5                 1 TRBV20-1*01             805        117            3
    #>  6                 1 TRBV20-1*01             379        117            3
    #>  7                 1 TRBV18*01               262        262            3
    #>  8                 1 TRBV18*01               586        262            3
    #>  9                 2 TRBV29-1*01            1100        695            3
    #> 10                 2 TRBV29-1*01             695        695            3
    #> # ℹ 47 more rows

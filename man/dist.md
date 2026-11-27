

# Immune receptor distances

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

A family of functions to calculate <strong>sequence distances between
immune receptors</strong>. These methods help you explore receptor
similarity, select thresholds for Immunoglobulin clone assignment, and
construct sparse receptor-similarity graphs.

<h4>Available functions</h4>

The following methods are available.

<strong>1) Hamming distances (<code>dist_hamm</code>).</strong> Compute
an upper-triangular sparse edge table of Hamming distances between
receptor IDs represented in the receptor annotations in
<code>idata</code>. Comparisons are restricted to equal sequence lengths
and equal values of <code>by</code>. The data are not grouped by
repertoire automatically. For paired-chain receptor schemas, the
canonical locus column is added to <code>by</code> automatically, so
distances are computed separately per locus. Add it explicitly if you
have a non-AIRR-C standard locus column name.

## Usage

<pre><code class='language-R'>dist_hamm(
  idata,
  seq_col = "cdr3_aa",
  by = NULL,
  max_dist = NULL,
  min_sim = NULL,
  sample_n = 0L,
  sample_by = "&lt;auto&gt;",
  mode = c("pairs", "graph")
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>An immundata::ImmunData object.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="seq_col">seq_col</code></td>
<td>A single character string naming a sequence feature in the receptor
schema.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>A character vector or <code>NULL</code>. These optional annotation
columns must be equal within a receptor pair. For B-cell clone
assignment, a typical choice is <code>c("v_call", "j_call")</code>; add
<code>imd_repertoire_id</code> when comparisons must stay within
repertoires, or <code>subject_id</code> if your data has such a
dedicated column for subject-specific metadata. Locus grouping is also
always applied for paired-chain receptor schemas.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="max_dist">max_dist</code></td>
<td>A single nonnegative number, a keyed threshold table, or
<code>NULL</code>. Scalar values in <code
style="white-space: pre;">(0, 1)</code> bound normalized Hamming
distance; integer values which are <code
style="white-space: pre;">&gt;=1</code> bound raw Hamming distance.
<code>NULL</code> (the default) and the legacy scalar value
<code>0</code> disable bounding. A threshold table must contain one or
more key columns plus a numeric <code>threshold</code> column strictly
inside <code style="white-space: pre;">(0, 1)</code>. Its key columns
must be included in <code>by</code> and must exactly cover the
corresponding input key combinations. Mutually exclusive with
<code>min_sim</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="min_sim">min_sim</code></td>
<td>A single number in <code style="white-space: pre;">[0, 1]</code> or
<code>NULL</code>, giving the minimum normalized Hamming similarity. It
is mutually exclusive with <code>max_dist</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="sample_n">sample_n</code></td>
<td>A nonnegative whole number giving the maximum number of distinct
receptor IDs to sample before distance calculation within each
<code>sample_by</code> group. <code>0</code> (the default) disables
sampling.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="sample_by">sample_by</code></td>
<td>A single character string or <code>NULL</code> selecting the column
to sample by. The default, <code>"&lt;auto&gt;"</code>, uses the first
available column among <code>subject_id</code>, <code>donor</code>,
<code>donor_id</code>, <code>patient_id</code>, <code>patient</code>,
<code>subject</code>, and <code>imd_repertoire_id</code>. The resolved
column must also be supplied through <code>by</code>. Set to
<code>NULL</code> to sample from the whole dataset without accounting
for groups. This argument is ignored when
<code>sample_n = 0</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="mode">mode</code></td>
<td>Output representation. <code>"pairs"</code> (the default) returns
the receptor-pair table described below. <code>"graph"</code> returns
contracted sequence-node edges and receptor-to-node membership. Grouping
columns cannot use names reserved for structural columns in the selected
output.</td>
</tr>
</tbody>
</table>

## Value

<h4>1) Hamming distances (<code>dist_hamm</code>)</h4>

With <code>mode = "pairs"</code>, a lazy duckplyr table with:

<ul>
<li> grouping columns supplied through <code>by</code>

</li>
<li> <code>locus</code> for paired-chain receptor schemas

</li>
<li> <code>imd_receptor_id_1</code>, <code>imd_receptor_id_2</code> &ndash; canonical receptor identifiers,
with <code>imd_receptor_id_1 &lt; imd_receptor_id_2</code>

</li>
<li> <code>seq_len</code>

</li>
<li> <code>dist</code> &ndash; raw Hamming distance

</li>
<li> <code>norm_dist</code> &ndash; <code>dist / seq_len</code>

</li>
<li> <code>sim</code> &ndash; normalized Hamming similarity, <code>1 - norm_dist</code>

</li></ul>

With <code>mode = "graph"</code>, an
<code>immunarch_sequence_graph</code> object containing two lazy
duckplyr tables. <code>edges</code> contains the effective grouping
columns, globally unique <code>from</code> and <code>to</code> node IDs,
and <code>seq_len</code>, <code>dist</code>, <code>norm_dist</code>, and
<code>sim</code>. <code>membership</code> contains the effective
grouping columns, the canonical receptor ID, and <code>node_id</code>. A
receptor can have multiple membership rows, and membership includes
nodes without edges. Nodes represent distinct combinations of the full
effective comparison group and a non-missing, non-empty sequence;
sequences themselves are not returned.

Call <code>dplyr::collect()</code> only when a returned table is small
enough for memory.

## Visualisation

Distance results can be passed directly to <code>vis()</code>.

<h4>1) Hamming distances (<code>dist_hamm</code>)</h4>

<code>vis()</code> plots the normalized distance to the nearest
non-identical neighbor for each receptor by default. This is the
distribution used for manual IG clone-threshold selection. When a
subject column can be inferred using the same priority as
<code>sample_by = "&lt;auto&gt;"</code>, one panel is drawn per subject;
otherwise, a single pooled distribution is drawn.

Use <code>mode = "all"</code> to plot every pairwise distance instead of
nearest neighbors. Pass one of the following values to <code>xval</code>
to change what will be plotted: <code>"norm_dist"</code> (the default),
<code>"dist"</code> (non-normalised distnaces), or <code>"sim"</code>
(similarity). Set <code>facet = NULL</code> to have one plot per all
input samples, or name a result column to facet by it.
<code>binwidth</code> and <code>title</code> customize the histogram.

## See Also

immundata::ImmunData

## Examples

``` r
library("immunarch")


idata <- get_test_idata()


# Build the computing query
res <- dist_hamm(idata)
# Check how results would look like
res
# Run the computation query and upload the results to the memory
res <- res |> collect()

# Bound Hamming distance.
# Non-normalised -- no receptors with more than 1 mismatch in the output
res <- dist_hamm(idata, max_dist = 1) |> collect()
res


# Normalised -- no receptors with more than .2 normalised distance in the output
res <- dist_hamm(idata, max_dist = .2) |> collect()
res

# Similarity, inverse to normalised distance -- no receptors with
# less than 0.9 similarity in the output
# The results are the same as for the previous normalised run
res <- dist_hamm(idata, min_sim = .8) |> collect()
res

# Plot normalized nearest-non-identical-neighbor distances.
vis(dist_hamm(idata))
# Plot all raw pairwise distances instead.
vis(dist_hamm(idata), mode = "all", xval = "dist")

#
# Typical use case for BCR data
#

# BCR-style V/J/length blocking with at most three mismatches, inside patients
clone_edges <- dist_hamm(
  idata_bcr,
  seq_col = "junction",
  by = c("patient", "v_call", "j_call"),
  max_dist = 3
)
```

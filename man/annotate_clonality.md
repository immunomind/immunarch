

# Annotate clonality - per-receptor labels for overabundance

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

A small family of helpers that <strong>add clonality labels to each
receptor</strong> in an immundata::ImmunData object.

<h4>Available functions</h4>
<ul>
<li> <code>annotate_clonality_rank()</code> - label by <strong>rank bins</strong> within each repertoire.

</li>
<li> <code>annotate_clonality_prop()</code> - label by <strong>proportion bins</strong> (named thresholds).

</li></ul>

<code>annotate_clonality_rank()</code> - for each repertoire, receptors
are ordered by within-repertoire abundance (proportion) and assigned a
<strong>rank bin</strong> label.

<code>annotate_clonality_prop()</code> - label each receptor by
<strong>proportion bin</strong> using named thresholds (matched in
descending order; else <code>"Ultra-rare"</code>).

## Usage

<pre><code class='language-R'>annotate_clonality_rank(idata, bins = c(10, 30, 100, 300, 1000, 10000, 1e+05))

annotate_clonality_prop(
  idata,
  bins = c(Hyperexpanded = 0.01, Large = 0.001, Medium = 1e-04, Small = 1e-05, Rare =
    1e-06)
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
id="bins">bins</code></td>
<td>A named numeric vector of thresholds (e.g.,
<code>c(Hyperexpanded = 1e-2, Large = 1e-3, ...)</code>). Names become
bin labels and must be non-empty. Internally sorted in descending order
for assignment. Result rows follow the supplied bin order, with
<code>"Ultra-rare"</code> last.</td>
</tr>
</tbody>
</table>

## Value

An immundata::ImmunData whose
<code style="white-space: pre;">$annotations</code> gains:

<ul>
<li> <code>clonal_rank_bin</code> - integer-like label with the applied rank threshold
(outside all thresholds -&gt; <code>NA</code>).

</li></ul>

An immundata::ImmunData whose
<code style="white-space: pre;">$annotations</code> gains:

<ul>
<li> <code>clonal_prop_bin</code> - label from <code>names(bins)</code> or <code>"Ultra-rare"</code>.

</li></ul>

## See Also

<ul>
<li> Per-repertoire summaries: airr_clonality

</li>
<li> Data container: immundata::ImmunData

</li></ul>

## Examples

``` r
library("immunarch")

idata <- get_test_idata() |> agg_repertoires("Therapy")
idata_rank <- annotate_clonality_rank(idata)
idata_prop <- annotate_clonality_prop(idata)
```

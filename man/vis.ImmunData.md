

# Visualise public receptor selection, incidence, or pairwise frequencies

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/v1_public_vis.R#L239)

## Description

Visualise selected receptor coverage and abundance with <code>plot =
"public"</code>, or receptor frequencies shared between two repertoires
with the default <code>plot = "pairwise"</code>. For
<code>annotate_public()</code> results, <code>plot = "incidence"</code>
shows receptor counts by the number of repertoires in which they occur.

In pairwise plots, receptor proportions are shown on log-scaled axes.
Point fill shows the log proportion ratio between the x and y
repertoires.

## Usage

<pre><code class='language-R'>## S3 method for class 'ImmunData'
vis(
  .data,
  x = 1,
  y = 2,
  density = TRUE,
  lm = TRUE,
  plot = c("pairwise", "public", "incidence"),
  reference = NULL,
  xval = "imd_repertoire_id",
  fill = NULL,
  title = "Selected immune receptors across repertoires",
  ...
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id=".data">.data</code></td>
<td>An immundata::ImmunData object with aggregated repertoires.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="x">x</code></td>
<td>Repertoire for the x-axis. Either a numeric repertoire index,
repertoire ID, or a repertoire schema value.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="y">y</code></td>
<td>Repertoire for the y-axis. Either a numeric repertoire index,
repertoire ID, or a repertoire schema value.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="density">density</code></td>
<td>Logical. If <code>TRUE</code>, show marginal distributions of unique
and public receptor proportions.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="lm">lm</code></td>
<td>Logical. If <code>TRUE</code>, add a linear model fit and adjusted
R-squared.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="plot">plot</code></td>
<td>One of <code>"pairwise"</code> (the default), <code>"public"</code>,
or <code>"incidence"</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="reference">reference</code></td>
<td>For <code>plot = "public"</code>, the original
<code>ImmunData</code> supplied to <code>vdj_public()</code>. It defines
the complete repertoire axis, including repertoires without selected
receptors. Required for this plot.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="xval">xval</code></td>
<td>For <code>plot = "public"</code>, the repertoire ID or metadata
column used for x-axis labels. Each value must identify one
repertoire.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="fill">fill</code></td>
<td>For <code>plot = "public"</code>, a repertoire metadata column used
for colour. <code>NULL</code> uses <code>strata_name</code> when strata
are defined; otherwise the plot uses a single colour.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="title">title</code></td>
<td>For <code>plot = "public"</code> or <code>"incidence"</code>, the
plot title, or <code>NULL</code> to omit it. Each view has its own
default title.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="...">...</code></td>
<td>Must be empty.</td>
</tr>
</tbody>
</table>

## Details

The incidence view stacks within-stratum receptor counts when strata are
defined. An identity observed in several strata contributes once to each
stratum's distribution; stack totals are not globally distinct receptor
counts. Without strata, one unstacked global distribution is shown and
every identity contributes once. This view uses statistics from the last
<code>annotate_public()</code> call; rerun annotation after changing the
population. Both publicness views use vertical legends on the right and
wrap long colour labels for display without changing the underlying
names.

## Value

A ggplot2 or patchwork object.

## Examples

``` r
library("immunarch")

idata |> vis()
idata |> vis(x = 1, y = 2)
idata |> vis(x = "Blood", y = "Tumor")
idata |> vis(x = 1, y = 2, density = FALSE)
selected <- vdj_public(idata)
vis(selected, plot = "public", reference = idata)
annotate_public(idata) |> vis(plot = "incidence")
```

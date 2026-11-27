

# Annotate receptors with publicness statistics

[**Source code**](https://github.com/immunomind/immunarch/tree/dev/R/v1_public_annotate.R#L61)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

## Usage

<pre><code class='language-R'>annotate_public(idata)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>An immundata::ImmunData object with aggregated repertoires.</td>
</tr>
</tbody>
</table>

## Details

Adds per-receptor publicness metrics to an immundata::ImmunData object.
Global metrics are always computed. If
<code style="white-space: pre;">$repertoires</code> contains
<code>imd_strata_id</code> (created by
<code>immundata::agg_strata()</code>), strata-specific metrics are added
as well. If strata is missing, the function reports that only global
metrics are produced.

For large datasets, this can be a computationally heavy operation. After
annotation, consider persisting the result as a snapshot with
<code>immundata::write_immundata()</code> to avoid recomputation.

## Value

An immundata::ImmunData with added annotation columns:

<ul>
<li> global metrics:
<code>imd_public_incidence</code>, <code>imd_public_incidence_prop</code>,
<code>imd_public_count_min</code>, <code>imd_public_count_max</code>,
<code>imd_public_count_mean</code>, <code>imd_public_count_median</code>,
<code>imd_public_prop_min</code>, <code>imd_public_prop_max</code>,
<code>imd_public_prop_mean</code>, <code>imd_public_prop_median</code>;

</li>
<li> if <code>imd_strata_id</code> exists in <code style="white-space: pre;">\$repertoires</code>, per-strata columns with suffix
<code style="white-space: pre;">_strata_{id}</code> for the same metric family.

</li></ul>

## Visualisation

Use <code>vis(annotated, plot = "incidence")</code> to show how many
receptor identities occur in exactly one, two, or more repertoires. With
strata, the plot stacks their within-stratum incidence distributions. A
receptor observed in several strata contributes once to each; the stack
heights sum within-stratum counts rather than globally distinct
identities. Without strata, the plot shows a single unstacked
distribution across all repertoires and omits the stratum legend. Long
stratum labels wrap in a vertical legend on the right; their underlying
names are unchanged. <code>title = NULL</code> omits the plot title.

Incidence annotations describe the population at the time of annotation.
Rerun <code>annotate_public()</code> after filtering or changing
repertoires or strata before plotting the new population.

## See Also

<code>vis()</code>, <code>vdj_public()</code>

## Examples

``` r
library("immunarch")

idata_global <- get_test_idata() |>
  agg_repertoires(c("Response", "Therapy")) |>
  annotate_public()
vis(idata_global, plot = "incidence")

idata_strata <- get_test_idata() |>
  agg_repertoires(c("Response", "Therapy")) |>
  agg_strata(schema = "Response") |>
  annotate_public()
vis(idata_strata, plot = "incidence")

write_immundata(idata_global, "snapshots/public-annotated-global")
write_immundata(idata_strata, "snapshots/public-annotated-strata")
```

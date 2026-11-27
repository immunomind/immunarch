

# One function to visualise them all

[**Source code**](https://github.com/immunomind/immunarch/tree/dev/R/v1_aaa_vis.R#L65)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

<strong><code>vis()</code> is a lightweight, quick-look plotting
helper.</strong> It's designed to help you <strong>visualise results
fast</strong> with sensible defaults. It automatically detects the input
type and chooses an appropriate visualisation (e.g., output from
<code>airr_desc_genes()</code> is recognised as gene usage values and
plotted without extra arguments).

<strong>vis() is not intended for publication-quality figures.</strong>
For serious, highly customised, or publication-ready plots, I recommend
building your graphics directly with <strong>ggplot2</strong>.

Information on how visualise a specific result from a domain, e.g.,
clonality or diversity, is available at each domain's help page in the
visualisation section. Take a look at the "See also" section down below
to <em>see</em> links to the domains, which support visualisations.

## Usage

<pre><code class='language-R'>vis(.data, ...)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id=".data">.data</code></td>
<td>The output from any immunarch analysis function. The function
automatically resolves to a correct visualisation.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="...">...</code></td>
<td>Any other arguments, see the "Details" section for specific
visualisation functions.</td>
</tr>
</tbody>
</table>

## Details

<code>vis()</code> chooses a plot based on the analysis result. Start
with <code>vis(result)</code>, then consult the documentation for the
analysis that created it to see its plotting defaults and supported
arguments.

E.g., for clonality results, see airr_clonality by running
<code>?airr_clonality</code>. Or see repsim For repertoire similarity by
running <code>?repsim</code>.

## Value

A ggplot2 object.

## See Also

<ul>
<li> airr_desc for descriptive statistics visualisations.

</li>
<li> airr_clonality for clonality visualisations.

</li>
<li> airr_diversity for diversity visualisations.

</li>
<li> repsim for repertoire similarity visualisations.

</li>
<li> dist for receptor distance visualisations.

</li>
<li> vdj_track for receptor tracking visualisations with <code>output = "df"</code>.

</li>
<li> vdj_public for selected receptor coverage and abundance. Pass its
<code>output = "df"</code> result and the original <code>reference</code> to <code>vis()</code>, or use
<code>plot = "public"</code> for its <code>ImmunData</code> result.

</li>
<li> annotate_public for incidence distributions with <code>plot = "incidence"</code>,
stacked by stratum when strata are defined and unstacked otherwise.

</li></ul>

## Examples

``` r
library("immunarch")

immdata <- get_test_idata() |> agg_repertoires("Therapy")
airr_desc_genes(immdata, gene_col = "v_call") |> vis()
repsim_jaccard(immdata) |> vis()
airr_diversity_pielou(immdata) |> vis()
airr_diversity_chao1(immdata) |> vis()
# Tracking summaries include explicit zeros for unobserved combinations.
vdj_track(immdata, n = 5, output = "df") |> vis()
vdj_public(immdata, output = "df") |> vis(reference = immdata)
annotate_public(immdata) |> vis(plot = "incidence")
```

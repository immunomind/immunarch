

# Describe one publicness selection rule

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/v1_public_vdj.R#L248)

## Description

Define a comparison population and inclusive incidence and abundance
bounds for <code>vdj_public()</code>. The constructor stores a rule
without evaluating data.

## Usage

<pre><code class='language-R'>make_public_filter(
  repertoires = NULL,
  strata = NULL,
  incidence_by = c("repertoire", "strata"),
  min_incidence = NULL,
  max_incidence = NULL,
  min_count = NULL,
  max_count = NULL,
  min_proportion = NULL,
  max_proportion = NULL
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="repertoires">repertoires</code></td>
<td>A numeric vector of repertoire IDs, or <code>NULL</code> for
all.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="strata">strata</code></td>
<td>Numeric strata IDs or existing <code>strata_name</code> labels, or
<code>NULL</code> for all. When both selectors are supplied, use their
intersection.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="incidence_by">incidence_by</code></td>
<td>A character scalar: <code>"repertoire"</code> counts distinct
repertoires; <code>"strata"</code> counts existing strata.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="min_incidence">min_incidence</code>, <code
id="max_incidence">max_incidence</code></td>
<td>Non-negative whole-number bounds, or <code>NULL</code>. Without
either bound, the minimum is one. An explicit maximum alone imposes no
minimum, so zero incidence can pass.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="min_count">min_count</code>, <code
id="max_count">max_count</code></td>
<td>Numeric count bounds, or <code>NULL</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="min_proportion">min_proportion</code>, <code
id="max_proportion">max_proportion</code></td>
<td>Numeric within-repertoire proportion bounds, or
<code>NULL</code>.</td>
</tr>
</tbody>
</table>

## Details

Abundance bounds qualify individual receptor occurrences before
incidence is counted. Repeated chain/cell abundance values are
deduplicated, never summed. For example,
<code style="white-space: pre;">max_count = 10, min_incidence = 2</code>
requires qualifying occurrences in two repertoires, allowing larger
occurrences elsewhere. <code style="white-space: pre;">min_count = 5,
max_incidence = 0</code> allows observed counts below five. Unknown
selectors and empty comparison populations are errors during evaluation.
Selectors do not match arbitrary metadata values.

## Value

A list describing one rule for <code>vdj_public()</code>.

## See Also

<code>vdj_public()</code> for strata preparation and applying several
rules.

## Examples

``` r
library("immunarch")

# Count >= 5 in at least one repertoire.
abundant <- make_public_filter(min_count = 5)


idata <- get_test_idata() |>
  agg_repertoires(c("Response", "Therapy"))

vdj_public(idata, abundant, output = "df") |>
  dplyr::select(cdr3_aa, v_call, Therapy, imd_orig_count) |>
  head()
```

    #> # A tibble: 6 × 4
    #>   cdr3_aa          v_call      Therapy imd_orig_count
    #>   <chr>            <chr>       <chr>            <dbl>
    #> 1 ASSPGQGLDTQY     TRBV18*01   ICI                 12
    #> 2 SVPGYYNEQF       TRBV29-1*01 ICI                  6
    #> 3 ASTVPGQGSLHTDTQY TRBV10-2*01 ICI                  8
    #> 4 ASSPGQGLDTQY     TRBV18*01   ICI                 12
    #> 5 SVLPEGGELF       TRBV29-1*01 ICI                  5
    #> 6 SAPRTDYGYT       TRBV20-1*01 ICI                  8

``` r
# At least one occurrence with count <= 10. Larger occurrences are allowed.
low_count <- make_public_filter(max_count = 10)

# No occurrence with count >= 5. Smaller counts are allowed.
below_threshold <- make_public_filter(min_count = 5, max_incidence = 0L)
```

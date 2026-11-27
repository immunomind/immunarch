

# Annotate immunoglobulin somatic hypermutation

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/#L)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

Calculate overall, V-region, and AIRR-region somatic-hypermutation
counts and frequencies from aligned observed and reconstructed germline
sequences.

The input must contain the columns selected by <code>sequence_col</code>
and <code>germline_col</code>, <code>v_sequence_start</code>,
<code>v_alignment_start</code>, and the canonical AIRR
<code style="white-space: pre;">\*\_start</code> and
<code style="white-space: pre;">\*\_end</code> columns for FWR1, CDR1,
FWR2, CDR2, FWR3, CDR3, and FWR4. The default germline column is
normally created with <code>annotate_ig_germlines()</code>.

## Usage

<pre><code class='language-R'>annotate_ig_shm(
  idata,
  sequence_col = "sequence_alignment",
  germline_col = "germline_alignment",
  by = "&lt;auto&gt;"
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>immundata::ImmunData. IgBLAST AIRR data to annotate.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="sequence_col">sequence_col</code></td>
<td>Character scalar. Aligned observed-sequence column.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="germline_col">germline_col</code></td>
<td>Character scalar. Aligned reconstructed germline column.
<code>sequence_col</code> and <code>germline_col</code> may select
alternative alignment columns, but those columns must use the same
alignment and coordinate system as the canonical AIRR region-coordinate
columns.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>Character vector, <code>"&lt;auto&gt;"</code>, or <code>NULL</code>.
Columns used to calculate group-level <code>mu_freq</code> summaries.
The default uses <code>clone_id</code>, with <code>subject_id</code>
when available; <code>locus</code> is added when available.
<code>NULL</code> skips grouped summaries.</td>
</tr>
</tbody>
</table>

## Details

AIRR region coordinates are interpreted as 1-based closed intervals on
the unaligned query sequence. <code>v_alignment_start</code> identifies
the aligned position corresponding to <code>v_sequence_start</code>.
While mapping query positions, <code>N</code> consumes a query
coordinate, whereas <code>.</code> and <code>-</code> do not.

A position is comparable only when both aligned characters are
unambiguous A, C, G, or T nucleotides. Thus <code>N</code>,
<code>.</code>, and <code>-</code> positions contribute to neither
mutation counts nor frequency denominators. <code>mu_count</code> is the
sum of classifiable replacement and silent mutations across all
available regions; <code>mu_freq</code> divides that count by the number
of comparable positions in those regions. Overall and combined CDR/FWR
replacement and silent frequencies use the sum of comparable positions
across their selected regions, not aggregated regional frequencies.
Individual-region frequencies use the comparable positions in that
region as their denominator.

Replacement and silent consequences are classified independently for
each mutation against the germline codon. Stop-producing substitutions
are not counted as either class. If one boundary of a region is missing,
all four outputs for that region are missing while other regions and
overall outputs remain available. A valid region with no comparable
positions has zero counts and missing frequencies. Combined CDR and
framework metrics use the available regions in their respective sets. If
no region in a set has both boundaries, all combined outputs for that
set are missing.

When grouping is enabled, <code>mu_freq_median</code> and
<code>mu_freq_mean</code> summarize non-missing <code>mu_freq</code>
values without weighting, and <code>mu_freq_n</code> counts their
contributors. Rows missing any grouping key receive missing median and
mean values and a zero count. Explicit <code>by</code> columns are used
exactly as supplied; include <code>locus</code> explicitly to keep heavy
and light chains separate. Use <code>dplyr::distinct()</code> on the
grouping columns and summaries when a group-level rather than row-level
result is needed.

## Value

immundata::ImmunData. The input rows with combined <code>mu_count</code>
and <code>mu_freq</code> columns plus replacement and silent counts and
frequencies overall, across the V region (FWR1 through FWR3), across all
CDRs, across all framework regions, and for each individual region. When
grouping is enabled, the result also includes
<code>mu_freq_median</code>, <code>mu_freq_mean</code>, and
<code>mu_freq_n</code> on every row.

## Why SHM results may differ between tools

Results can differ between SHM tools, even for the same sequences. This
function counts a mutation only when its codon (a group of three bases)
is complete and has no alignment gap. Other tools may count mutations in
codons with gaps. Removing <code>.</code> characters from an IMGT
alignment can also shift codon boundaries and change which mutations are
counted, including mutations that create a stop codon. For a fair
comparison, use matching alignments and region positions, and check how
each tool handles gaps and stop codons.

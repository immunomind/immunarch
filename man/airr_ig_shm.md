

# Summarize immunoglobulin somatic hypermutation per repertoire

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/#L)

## Description

Summarize the mutation counts and frequencies produced by
<code>annotate_ig_shm()</code>. Aggregation runs in DuckDB on the
duckplyr-backed input.

## Usage

<pre><code class='language-R'>airr_ig_shm(
  idata,
  by = NA_character_,
  clone_col = "clone_id",
  autojoin = getOption(IMMUNARCH_AUTOJOIN_OPTION, IMMUNARCH_AUTOJOIN_DEFAULT),
  format = c("long", "wide")
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>An immundata::ImmunData object with SHM annotations and clones.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>Character vector, <code>NA</code>, or <code>NULL</code>. The default
uses the canonical locus column when available. Use <code>NULL</code> to
pool loci. Repertoires are always summarized separately; clone IDs
cannot be included in <code>by</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="clone_col">clone_col</code></td>
<td>Column containing clone assignments. Missing clone IDs are included
in sequence summaries but excluded from clone summaries.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="autojoin">autojoin</code></td>
<td>Logical. If TRUE, join repertoire metadata by the schema repertoire
id. Change the default behaviour by calling
<code>options(immunarch.autojoin = FALSE)</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="format">format</code></td>
<td>Inherited from the common AIRR interface. This method always returns
long-format summaries.</td>
</tr>
</tbody>
</table>

## Details

This function uses the mutation counts and frequencies already
calculated by <code>annotate_ig_shm()</code>. Here, a sequence means one
immunoglobulin chain sequence in the input.

The function summarizes these counts and frequencies in two ways:

<ul>
<li> Across sequences: each sequence has equal weight. A clone with more
sequences therefore contributes more to the result.

</li>
<li> Across clones: first calculate a mean or median within each clone, then
summarize these clone means or medians. Each clone has equal weight,
regardless of how many sequences it contains.

</li></ul>

These calculations do not use read counts or UMI counts as weights.

## Value

A long-format table of repertoire statistics. Each row contains the
repertoire ID, any <code>by</code> columns, <code>region</code>,
<code>mutation</code>, <code>measure</code>, and the 24 statistics
listed above. Repertoire metadata is added by default.

## Output statistics

Each repertoire has several output rows. The columns
<code>region</code>, <code>mutation</code>, and <code>measure</code>
tell you which mutation counts or frequencies the row describes:

<table data-quarto-postprocess="true">
<tbody>
<tr>
<td style="text-align: left;">Column</td>
<td style="text-align: left;">Meaning</td>
</tr>
<tr>
<td style="text-align: left;"><code>region</code></td>
<td style="text-align: left;">The region being analyzed:
<code>overall</code> for all available annotated regions; <code>v</code>
for FWR1 through FWR3; <code>cdr</code> for all available CDRs;
<code>fwr</code> for all available framework regions; or one region such
as <code>cdr3</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>mutation</code></td>
<td style="text-align: left;"><code>replacement</code> for mutations
that change the amino acid, <code>silent</code> for mutations that do
not change the amino acid, or <code>total</code> for replacement plus
silent mutations. See <code>annotate_ig_shm()</code> for counting
rules.</td>
</tr>
<tr>
<td style="text-align: left;"><code>measure</code></td>
<td style="text-align: left;"><code>count</code> for the number of
mutations, or <code>freq</code> for mutation frequency. Frequencies are
fractions: <code>0.05</code> means 5%.</td>
</tr>
</tbody>
</table>

The other columns contain the statistics for that row. For example,
consider four sequences in one repertoire. Clone A has total CDR3
mutation counts of 0, 0, and 9. Clone B has one sequence with 6 total
CDR3 mutations. Part of the output would be:

<table data-quarto-postprocess="true">
<tbody>
<tr>
<td style="text-align: left;">imd_repertoire_id</td>
<td style="text-align: left;">region</td>
<td style="text-align: left;">mutation</td>
<td style="text-align: left;">measure</td>
<td style="text-align: right;">sequence_mean</td>
<td style="text-align: right;">clone_mean_mean</td>
<td style="text-align: right;">clone_median_mean</td>
</tr>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;">cdr3</td>
<td style="text-align: left;">total</td>
<td style="text-align: left;">count</td>
<td style="text-align: right;">3.75</td>
<td style="text-align: right;">4.5</td>
<td style="text-align: right;">3</td>
</tr>
</tbody>
</table>

All three statistics in this row describe total CDR3 mutation counts:

<ul>
<li> <code>sequence_mean</code>: average the four counts: <code>(0 + 0 + 9 + 6) / 4 = 3.75</code>.

</li>
<li> <code>clone_mean_mean</code>: the clone means are 3 and 6. Their mean is 4.5.

</li>
<li> <code>clone_median_mean</code>: the clone medians are 0 and 6. Their mean is 3.

</li></ul>

A row with <code>measure = "freq"</code> contains the same kinds of
statistics for mutation frequencies instead of mutation counts. Results
are calculated separately for each repertoire and each group selected
with <code>by</code>, such as locus. Repertoire metadata is added when
<code>autojoin = TRUE</code>.

The full list of 24 statistics is below. In every description, the
mutation count or frequency is the one specified by <code>region</code>,
<code>mutation</code>, and <code>measure</code> in that row.
<code>q25</code> and <code>q75</code> are the lower and upper quartiles
(25th and 75th percentiles).

<table data-quarto-postprocess="true">
<tbody>
<tr>
<td style="text-align: left;">Column</td>
<td style="text-align: left;">Meaning</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_sequences</code></td>
<td style="text-align: left;">Number of input sequences, including those
with a missing mutation count or frequency, or a missing clone ID.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_valid</code></td>
<td style="text-align: left;">Number of sequences with an available
mutation count or frequency.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_missing</code></td>
<td style="text-align: left;">Number of sequences whose mutation count
or frequency is <code>NA</code>:
<code>n_sequences - n_valid</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_valid_without_clone</code></td>
<td style="text-align: left;">Number of sequences with an available
mutation count or frequency but no clone ID.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_zero</code></td>
<td style="text-align: left;">Number of sequences with a mutation count
or frequency of zero.</td>
</tr>
<tr>
<td style="text-align: left;"><code>zero_prop</code></td>
<td style="text-align: left;">Proportion of sequences with zero
mutations of this type in this region, among sequences with an available
count or frequency: <code>n_zero / n_valid</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>sequence_mean</code></td>
<td style="text-align: left;">Mean mutation count or frequency across
sequences. Each sequence has equal weight.</td>
</tr>
<tr>
<td style="text-align: left;"><code>sequence_median</code></td>
<td style="text-align: left;">Median mutation count or frequency across
sequences.</td>
</tr>
<tr>
<td style="text-align: left;"><code>sequence_q25</code></td>
<td style="text-align: left;">Lower quartile of mutation counts or
frequencies across sequences.</td>
</tr>
<tr>
<td style="text-align: left;"><code>sequence_q75</code></td>
<td style="text-align: left;">Upper quartile of mutation counts or
frequencies across sequences.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_clones</code></td>
<td style="text-align: left;">Number of clones with a known clone ID,
including clones for which all mutation counts or frequencies are
missing.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_valid_clones</code></td>
<td style="text-align: left;">Number of clones with at least one
sequence whose mutation count or frequency is available.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_valid_with_clone</code></td>
<td style="text-align: left;">Number of sequences with both an available
mutation count or frequency and a clone ID.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_mean_weighted_mean</code></td>
<td style="text-align: left;">Mean mutation count or frequency across
sequences with a clone ID. This is also the mean of clone means weighted
by the number of sequences with available counts or frequencies in each
clone.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_mean_mean</code></td>
<td style="text-align: left;">Calculate the mean mutation count or
frequency in each clone, then take the mean of these clone means. Each
clone has equal weight.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_mean_median</code></td>
<td style="text-align: left;">Calculate the mean mutation count or
frequency in each clone, then take the median of these clone means.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_mean_q25</code></td>
<td style="text-align: left;">Lower quartile of the clone means.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_mean_q75</code></td>
<td style="text-align: left;">Upper quartile of the clone means.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_median_mean</code></td>
<td style="text-align: left;">Calculate the median mutation count or
frequency in each clone, then take the mean of these clone medians. Each
clone has equal weight.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_median_median</code></td>
<td style="text-align: left;">Calculate the median mutation count or
frequency in each clone, then take the median of these clone
medians.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_median_q25</code></td>
<td style="text-align: left;">Lower quartile of the clone medians.</td>
</tr>
<tr>
<td style="text-align: left;"><code>clone_median_q75</code></td>
<td style="text-align: left;">Upper quartile of the clone medians.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_zero_clone_means</code></td>
<td style="text-align: left;">Number of clones whose mean mutation count
or frequency is zero.</td>
</tr>
<tr>
<td style="text-align: left;"><code>n_zero_clone_medians</code></td>
<td style="text-align: left;">Number of clones whose median mutation
count or frequency is zero.</td>
</tr>
</tbody>
</table>

Missing mutation counts or frequencies (<code>NA</code>) are excluded
from means, medians, quartiles, and the denominator of
<code>zero_prop</code>. Zero is included. Sequences without clone IDs
are included in sequence statistics but excluded from clone statistics.
Clones with no available mutation counts or frequencies are included in
<code>n_clones</code> but excluded from the clone means, medians, and
quartiles. If none of the required counts or frequencies are available,
the mean, median, quartiles, and <code>zero_prop</code> are
<code>NA</code>. Quartiles use linear interpolation.

Combined region frequencies come from <code>annotate_ig_shm()</code>,
which divides the total mutation count by the total number of comparable
bases in those regions. This function does not average the frequencies
of different regions to obtain a combined frequency. Regions without the
required input columns are skipped.

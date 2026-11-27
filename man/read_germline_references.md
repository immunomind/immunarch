

# Read IMGT immunoglobulin germline references

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/v1_ig_germlines_io.R#L26)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

V and J genes provide parts of an antibody sequence. To measure
mutations in a B-cell clone, we need to compare its sequences with
reference versions of those genes before somatic hypermutation.

This function reads nucleotide V and J genes from an IMGT FASTA file for
use by <code>annotate_ig_germlines()</code>. It keeps the IMGT gaps
(<code>.</code>) in the reference sequences.

## Usage

<pre><code class='language-R'>read_germline_references(path)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="path">path</code></td>
<td>Character scalar. Path to an IMGT FASTA file.</td>
</tr>
</tbody>
</table>

## Value

A tibble with <code>species</code>, <code>locus</code>,
<code>segment</code>, <code>allele_id</code>, and <code>sequence</code>
columns.

## Choosing files

Use IMGT-gapped V and J FASTA files from the same IMGT release used for
your IgBlast run. Read each needed file and combine the resulting
tables. Do not use the ungapped FASTA from the searchable IgBlast
database: this function does not add gaps or convert positions.

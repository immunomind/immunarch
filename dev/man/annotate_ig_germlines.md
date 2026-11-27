

# Reconstruct clone-level immunoglobulin germlines

[**Source code**](https://github.com/immunomind/immunarch/tree/dev/R/#L)

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

B cells in the same clone come from one starting cell, but their
antibody sequences can gain different mutations over time. To study
these mutations, we need a shared sequence to compare them with.

The method annotates each clone and locus with a representative germline
sequence from <code>references</code> – usually IMGT germline table.

The result can be used as the reference for somatic hypermutation
analysis with <code>annotate_ig_shm()</code>.

## Usage

<pre><code class='language-R'>annotate_ig_germlines(
  idata,
  references,
  species = NULL,
  by = "&lt;auto&gt;",
  sequence_id_col = "sequence_id",
  sequence_col = "sequence_alignment",
  v_col = "v_call",
  j_col = "j_call",
  v_start_col = "v_germline_start",
  v_end_col = "v_germline_end",
  j_start_col = "j_germline_start",
  j_end_col = "j_germline_end",
  germline_col = "germline_alignment"
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="idata">idata</code></td>
<td>immundata::ImmunData. AIRR data to annotate.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="references">references</code></td>
<td>A table returned by <code>read_germline_references()</code>
containing IMGT-gapped reconstruction references. These are not the
ungapped reference files used to build an IgBLAST database; see
<strong>Assumptions and reference compatibility</strong>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="species">species</code></td>
<td>Character scalar, or <code>NULL</code> to use the sole reference
species.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>Character vector or <code>"&lt;auto&gt;"</code>. Clone grouping
columns. Automatic grouping uses <code>clone_id</code>, with
<code>subject_id</code> when available. Locus is always added.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="sequence_id_col">sequence_id_col</code>, <code
id="sequence_col">sequence_col</code>, <code id="v_col">v_col</code>,
<code id="j_col">j_col</code></td>
<td>Character scalars selecting the sequence ID, IMGT-gapped aligned
sequence, V-call, and J-call columns.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="v_start_col">v_start_col</code>, <code
id="v_end_col">v_end_col</code>, <code
id="j_start_col">j_start_col</code>, <code
id="j_end_col">j_end_col</code></td>
<td>Character scalars selecting the closed, IMGT-gapped
germline-coordinate columns.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="germline_col">germline_col</code></td>
<td>Character scalar naming the reconstructed output column.</td>
</tr>
</tbody>
</table>

## Value

immundata::ImmunData with one clone/locus germline assigned to each
input row. Groups which cannot be reconstructed receive <code>NA</code>.

## Assumptions and reference compatibility

This function reconstructs germlines by indexing directly into the
supplied reference strings. It does not insert IMGT numbering gaps,
convert coordinates, or otherwise align the input sequences and
references. Consequently, all of the following assumptions must hold:

<ul>
<li> <code>sequence_col</code> contains IMGT-gapped alignments, including the standard <code>.</code>
numbering spacers.

</li>
<li> The germline start and end columns use the same IMGT-gapped coordinate
system as <code>sequence_col</code> and <code>references</code>.

</li>
<li> V references in <code>references</code> contain IMGT numbering gaps. In particular,
do not supply the ungapped FASTA files used to build an IgBLAST database.

</li>
<li> The reference alleles and release match those used during assignment of the
V and J calls. Missing or different alleles cannot be reconstructed safely.

</li></ul>

An Immcantation-style workflow normally maintains two representations of
the same reference release: ungapped sequences for the IgBLAST database
and IMGT-gapped sequences for germline reconstruction. This function
requires the latter. Matching references therefore means the same allele
set and release, not the literal ungapped FASTA passed to IgBLAST.

<code>read_germline_references()</code> preserves gaps that are already
present in a FASTA file, but neither it nor this function converts
ungapped references to IMGT-gapped references. Supplying incompatible
ungapped data can produce a plausible-length but incorrectly shifted
germline without an error.

The reconstructed value corresponds conceptually to Dowser's
<code>germline_alignment_d_mask</code>: the retained V and J portions
surround an <code>N</code>-masked middle region. This function does not
reconstruct a full V(D)J germline or infer a D segment.

## Preparing the input

To create the BCR repertoire, first you need to run IgBlast either
standalone or through <code>immundapter::tool_igblast()</code> with an
IMGT database on your FASTQ files or 10XGenomics contigs files.
<code>immundapter</code> automatically prepares the alignment and V
positions for IMGT-gapped references. Then, read V and J references from
the same IMGT release with <code>read_germline_references()</code>. See
the "Ig tutorial" on the website for the detailed workflow – it is a bit
complicated, as per usual with BCR analysis.

Other AIRR data can be used when its alignment and V/J positions match
the supplied IMGT-gapped references. But please mind that raw IgBlast
AIRR output has not been prepared this way. This function does not
convert positions or add gaps, so mismatched input can give an incorrect
sequence without an error.

## Output

This is a D-masked germline, not a complete ancestral receptor sequence.
By default, it replaces an existing <code>germline_alignment</code>
column. Use <code>germline_col</code> to keep that column.

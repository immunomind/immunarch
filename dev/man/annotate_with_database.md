

# Annotate ImmunData with an external sequence database (e.g. VDJdb)

[**Source code**](https://github.com/immunomind/immunarch/tree/dev/R/v1_database_annotate.R#L19)

## Description

Add annotations to ImmunData from an external database such as VDJdb.

## Usage

<pre><code class='language-R'>annotate_with_database(idata, db, by, label_col = "species", ...)
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
id="db">db</code></td>
<td>A data frame, tibble, <code>duckplyr_df</code>, or a single
character string containing a path to a delimited or Parquet file. Files
are read via duckplyr.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="by">by</code></td>
<td>A named character vector of join columns, of the form
<code>c(idata_col = "db_col")</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="label_col">label_col</code></td>
<td>A single character string naming the column from <code>db</code> to
add to <code>idata</code> (for VDJdb, this could be
<code>"species"</code>).</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="...">...</code></td>
<td>Optional dplyr or duckplyr filter expressions evaluated on
<code>db</code> before joining (e.g.
<code>species == "HomoSapiens"</code>, <code>gene == "TRB"</code>).</td>
</tr>
</tbody>
</table>

## Value

A new immundata::ImmunData object with the selected database label added
as an annotation column.

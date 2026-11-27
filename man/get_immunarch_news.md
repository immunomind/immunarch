

# Get the Latest immunarch Update

[**Source code**](https://github.com/immunomind/immunarch/tree/master/R/v1_migration_utils.R#L20)

## Description

Retrieves an update message for immunarch.

## Usage

<pre><code class='language-R'>get_immunarch_news(datepoint = "latest")
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="datepoint">datepoint</code></td>
<td>A single character string specifying the update date. Use
<code>"latest"</code> for the most recent update or supply a valid date
key (e.g., <code>"Apr 2025"</code>).</td>
</tr>
</tbody>
</table>

## Details

If <code>datepoint</code> is set to <code>"latest"</code>, the function
returns the most recent update. Otherwise, specify the update date key
(e.g., <code>"Apr 2025"</code>) to retrieve that particular update. If
no matching update is found, a warning is issued along with available
update keys.

## Value

A character string with the update details or a warning if the key is
not found.

## See Also

<code>list_immunarch_news()</code>

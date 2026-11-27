

# Analyse and annotate immunoglobulin receptors

## Description

<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="../help/figures/lifecycle-experimental.svg" alt='[Experimental]'></a>

A family of functions for analysing and adding biologically useful
immunoglobulin annotations to receptors. These annotations help connect
receptor sequence variation with clonal expansion and B-cell evolution.

Each function has its own dedicated help page, so click on the function
names to open corresponding pages.

<h4>Available functions</h4>
<ol>
<li> <strong>Clonal family inference (<code>annotate_ig_clones()</code>).</strong> Usually the very first step in
any of the Ig analysis workflows. The function clusters Ig receptors to clonal families using
heavy-chain or paired heavy-light junction similarity by Hamming distance.

</li>
<li> <strong>Germline annotation (<code>annotate_ig_germlines()</code>).</strong> The first step in the
somatic hypermutation (SHM) analysis workflows. This function annotates the
Ig receptors with information about germline genes for the subsequent SHM analysis.

</li>
<li> <strong>Somatic hypermutation computation (<code>annotate_ig_shm()</code>).</strong> The second and the final step in the
SHM analysis workflows. This function computes the SHM per-receptor and aggregated
statistics and annotates the receptors with them. After that you can use
<code>airr_ig_shm()</code> to explore and analyse the repertoire-wide and family-wide statistics
in a more convenient form.

</li></ol>

## See Also

<code>annotate_ig_clones()</code>, <code>annotate_ig_germlines()</code>,
<code>annotate_ig_shm()</code>

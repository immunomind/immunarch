#' Retired immunarch v0 functions
#'
#' These functions stopped working in immunarch 0.12.0. Calls signal a
#' deprecation error. Use immunarch 0.11.x to run an existing v0 workflow.
#' For migration steps, see the
#' [migration guide](https://immunomind.github.io/docs/tutorials/migration/).
#'
#' @param ... Arguments from a legacy call.
#' @return These functions signal a deprecation error.
#' @name immunarch-v0-defunct
#' @concept retired_api
NULL

stop_v0_api <- function(what, with = NULL) {
  lifecycle::deprecate_stop(
    "0.12.0",
    what,
    with = with,
    details = "See the migration guide at https://immunomind.github.io/docs/tutorials/migration/"
  )
}

#' @rdname immunarch-v0-defunct
#' @export
apply_asymm <- function(...) stop_v0_api("apply_asymm()")

#' @rdname immunarch-v0-defunct
#' @export
apply_symm <- function(...) stop_v0_api("apply_symm()")

#' @rdname immunarch-v0-defunct
#' @export
bunch_translate <- function(...) stop_v0_api("bunch_translate()")

#' @rdname immunarch-v0-defunct
#' @export
coding <- function(...) stop_v0_api("coding()")

#' @rdname immunarch-v0-defunct
#' @export
cross_entropy <- function(...) stop_v0_api("cross_entropy()")

#' @rdname immunarch-v0-defunct
#' @export
dbAnnotate <- function(...) stop_v0_api("dbAnnotate()", "annotate_with_database()")

#' @rdname immunarch-v0-defunct
#' @export
dbLoad <- function(...) stop_v0_api("dbLoad()")

#' @rdname immunarch-v0-defunct
#' @export
entropy <- function(...) stop_v0_api("entropy()")

#' @rdname immunarch-v0-defunct
#' @export
exclude <- function(...) stop_v0_api("exclude()")

#' @rdname immunarch-v0-defunct
#' @export
fixVis <- function(...) stop_v0_api("fixVis()")

#' @rdname immunarch-v0-defunct
#' @export
geneUsage <- function(...) stop_v0_api("geneUsage()", "airr_desc_genes()")

#' @rdname immunarch-v0-defunct
#' @export
geneUsageAnalysis <- function(...) stop_v0_api("geneUsageAnalysis()")

#' @rdname immunarch-v0-defunct
#' @export
gene_stats <- function(...) stop_v0_api("gene_stats()")

#' @rdname immunarch-v0-defunct
#' @export
getKmers <- function(...) stop_v0_api("getKmers()")

#' @rdname immunarch-v0-defunct
#' @export
get_genes <- function(...) stop_v0_api("get_genes()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_dbscan <- function(...) stop_v0_api("immunr_dbscan()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_hclust <- function(...) stop_v0_api("immunr_hclust()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_kmeans <- function(...) stop_v0_api("immunr_kmeans()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_mds <- function(...) stop_v0_api("immunr_mds()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_pca <- function(...) stop_v0_api("immunr_pca()")

#' @rdname immunarch-v0-defunct
#' @export
immunr_tsne <- function(...) stop_v0_api("immunr_tsne()")

#' @rdname immunarch-v0-defunct
#' @export
inc_overlap <- function(...) stop_v0_api("inc_overlap()")

#' @rdname immunarch-v0-defunct
#' @export
include <- function(...) stop_v0_api("include()")

#' @rdname immunarch-v0-defunct
#' @export
inframes <- function(...) stop_v0_api("inframes()")

#' @rdname immunarch-v0-defunct
#' @export
interval <- function(...) stop_v0_api("interval()")

#' @rdname immunarch-v0-defunct
#' @export
js_div <- function(...) stop_v0_api("js_div()")

#' @rdname immunarch-v0-defunct
#' @export
kl_div <- function(...) stop_v0_api("kl_div()")

#' @rdname immunarch-v0-defunct
#' @export
kmer_profile <- function(...) stop_v0_api("kmer_profile()")

#' @rdname immunarch-v0-defunct
#' @export
lessthan <- function(...) stop_v0_api("lessthan()")

#' @rdname immunarch-v0-defunct
#' @export
morethan <- function(...) stop_v0_api("morethan()")

#' @rdname immunarch-v0-defunct
#' @export
noncoding <- function(...) stop_v0_api("noncoding()")

#' @rdname immunarch-v0-defunct
#' @export
outofframes <- function(...) stop_v0_api("outofframes()")

#' @rdname immunarch-v0-defunct
#' @export
pubRep <- function(...) stop_v0_api("pubRep()")

#' @rdname immunarch-v0-defunct
#' @export
pubRepApply <- function(...) stop_v0_api("pubRepApply()")

#' @rdname immunarch-v0-defunct
#' @export
pubRepFilter <- function(...) stop_v0_api("pubRepFilter()")

#' @rdname immunarch-v0-defunct
#' @export
pubRepStatistics <- function(...) stop_v0_api("pubRepStatistics()")

#' @rdname immunarch-v0-defunct
#' @export
publicRepertoire <- function(...) stop_v0_api("publicRepertoire()")

#' @rdname immunarch-v0-defunct
#' @export
publicRepertoireApply <- function(...) stop_v0_api("publicRepertoireApply()")

#' @rdname immunarch-v0-defunct
#' @export
publicRepertoireFilter <- function(...) stop_v0_api("publicRepertoireFilter()")

#' @rdname immunarch-v0-defunct
#' @export
public_matrix <- function(...) stop_v0_api("public_matrix()")

#' @rdname immunarch-v0-defunct
#' @export
repAlignLineage <- function(...) stop_v0_api("repAlignLineage()")

#' @rdname immunarch-v0-defunct
#' @export
repClonalFamily <- function(...) stop_v0_api("repClonalFamily()")

#' @rdname immunarch-v0-defunct
#' @export
repClonality <- function(...) stop_v0_api("repClonality()")

#' @rdname immunarch-v0-defunct
#' @export
repDiversity <- function(...) stop_v0_api("repDiversity()")

#' @rdname immunarch-v0-defunct
#' @export
repExplore <- function(...) stop_v0_api("repExplore()")

#' @rdname immunarch-v0-defunct
#' @export
repFilter <- function(...) stop_v0_api("repFilter()", "immundata::filter_immundata()")

#' @rdname immunarch-v0-defunct
#' @export
repGermline <- function(...) stop_v0_api("repGermline()")

#' @rdname immunarch-v0-defunct
#' @export
repLoad <- function(...) stop_v0_api("repLoad()", "immundata::read_repertoires()")

#' @rdname immunarch-v0-defunct
#' @export
repOverlap <- function(...) stop_v0_api("repOverlap()")

#' @rdname immunarch-v0-defunct
#' @export
repOverlapAnalysis <- function(...) stop_v0_api("repOverlapAnalysis()")

#' @rdname immunarch-v0-defunct
#' @export
repSample <- function(...) stop_v0_api("repSample()")

#' @rdname immunarch-v0-defunct
#' @export
repSave <- function(...) stop_v0_api("repSave()", "immundata::write_immundata()")

#' @rdname immunarch-v0-defunct
#' @export
repSomaticHypermutation <- function(...) stop_v0_api("repSomaticHypermutation()")

#' @rdname immunarch-v0-defunct
#' @export
select_barcodes <- function(...) stop_v0_api("select_barcodes()")

#' @rdname immunarch-v0-defunct
#' @export
select_clusters <- function(...) stop_v0_api("select_clusters()")

#' @rdname immunarch-v0-defunct
#' @export
seqCluster <- function(...) stop_v0_api("seqCluster()")

#' @rdname immunarch-v0-defunct
#' @export
seqDist <- function(...) stop_v0_api("seqDist()")

#' @rdname immunarch-v0-defunct
#' @export
spectratype <- function(...) stop_v0_api("spectratype()")

#' @rdname immunarch-v0-defunct
#' @export
split_to_kmers <- function(...) stop_v0_api("split_to_kmers()")

#' @rdname immunarch-v0-defunct
#' @export
top <- function(...) stop_v0_api("top()")

#' @rdname immunarch-v0-defunct
#' @export
trackClonotypes <- function(...) stop_v0_api("trackClonotypes()", "vdj_track()")

#' @rdname immunarch-v0-defunct
#' @export
vis_bar <- function(...) stop_v0_api("vis_bar()")

#' @rdname immunarch-v0-defunct
#' @export
vis_box <- function(...) stop_v0_api("vis_box()")

#' @rdname immunarch-v0-defunct
#' @export
vis_circos <- function(...) stop_v0_api("vis_circos()")

#' @rdname immunarch-v0-defunct
#' @export
vis_heatmap <- function(...) stop_v0_api("vis_heatmap()")

#' @rdname immunarch-v0-defunct
#' @export
vis_heatmap2 <- function(...) stop_v0_api("vis_heatmap2()")

#' @rdname immunarch-v0-defunct
#' @export
vis_hist <- function(...) stop_v0_api("vis_hist()")

#' @rdname immunarch-v0-defunct
#' @export
vis_immunr_kmer_profile_main <- function(...) stop_v0_api("vis_immunr_kmer_profile_main()")

#' @rdname immunarch-v0-defunct
#' @export
vis_seqlogo <- function(...) stop_v0_api("vis_seqlogo()")

#' @rdname immunarch-v0-defunct
#' @export
vis_textlogo <- function(...) stop_v0_api("vis_textlogo()")

# Legacy result classes can still occur in saved analyses.

#' @export
vis.clonal_family <- function(.data, ...) stop_v0_api("vis.clonal_family()")

#' @export
vis.clonal_family_tree <- function(.data, ...) stop_v0_api("vis.clonal_family_tree()")

#' @export
vis.immunr_chao1 <- function(.data, ...) stop_v0_api("vis.immunr_chao1()")

#' @export
vis.immunr_clonal_prop <- function(.data, ...) stop_v0_api("vis.immunr_clonal_prop()")

#' @export
vis.immunr_dbscan <- function(.data, ...) stop_v0_api("vis.immunr_dbscan()")

#' @export
vis.immunr_div <- function(.data, ...) stop_v0_api("vis.immunr_div()")

#' @export
vis.immunr_dxx <- function(.data, ...) stop_v0_api("vis.immunr_dxx()")

#' @export
vis.immunr_dynamics <- function(.data, ...) stop_v0_api("vis.immunr_dynamics()")

#' @export
vis.immunr_exp_clones <- function(.data, ...) stop_v0_api("vis.immunr_exp_clones()")

#' @export
vis.immunr_exp_count <- function(.data, ...) stop_v0_api("vis.immunr_exp_count()")

#' @export
vis.immunr_exp_len <- function(.data, ...) stop_v0_api("vis.immunr_exp_len()")

#' @export
vis.immunr_exp_vol <- function(.data, ...) stop_v0_api("vis.immunr_exp_vol()")

#' @export
vis.immunr_gene_usage <- function(.data, ...) stop_v0_api("vis.immunr_gene_usage()")

#' @export
vis.immunr_ginisimp <- function(.data, ...) stop_v0_api("vis.immunr_ginisimp()")

#' @export
vis.immunr_gu_matrix <- function(.data, ...) stop_v0_api("vis.immunr_gu_matrix()")

#' @export
vis.immunr_hclust <- function(.data, ...) stop_v0_api("vis.immunr_hclust()")

#' @export
vis.immunr_hill <- function(.data, ...) stop_v0_api("vis.immunr_hill()")

#' @export
vis.immunr_homeo <- function(.data, ...) stop_v0_api("vis.immunr_homeo()")

#' @export
vis.immunr_inc_overlap <- function(.data, ...) stop_v0_api("vis.immunr_inc_overlap()")

#' @export
vis.immunr_invsimp <- function(.data, ...) stop_v0_api("vis.immunr_invsimp()")

#' @export
vis.immunr_kmeans <- function(.data, ...) stop_v0_api("vis.immunr_kmeans()")

#' @export
vis.immunr_kmer_profile_pfm <- function(.data, ...) stop_v0_api("vis.immunr_kmer_profile_pfm()")

#' @export
vis.immunr_kmer_profile_ppm <- function(.data, ...) stop_v0_api("vis.immunr_kmer_profile_ppm()")

#' @export
vis.immunr_kmer_profile_pwm <- function(.data, ...) stop_v0_api("vis.immunr_kmer_profile_pwm()")

#' @export
vis.immunr_kmer_profile_self <- function(.data, ...) stop_v0_api("vis.immunr_kmer_profile_self()")

#' @export
vis.immunr_kmer_table <- function(.data, ...) stop_v0_api("vis.immunr_kmer_table()")

#' @export
vis.immunr_mds <- function(.data, ...) stop_v0_api("vis.immunr_mds()")

#' @export
vis.immunr_ov_matrix <- function(.data, ...) stop_v0_api("vis.immunr_ov_matrix()")

#' @export
vis.immunr_pca <- function(.data, ...) stop_v0_api("vis.immunr_pca()")

#' @export
vis.immunr_public_repertoire <- function(.data, ...) stop_v0_api("vis.immunr_public_repertoire()")

#' @export
vis.immunr_public_statistics <- function(.data, ...) stop_v0_api("vis.immunr_public_statistics()")

#' @export
vis.immunr_rare_prop <- function(.data, ...) stop_v0_api("vis.immunr_rare_prop()")

#' @export
vis.immunr_rarefaction <- function(.data, ...) stop_v0_api("vis.immunr_rarefaction()")

#' @export
vis.immunr_spectr <- function(.data, ...) stop_v0_api("vis.immunr_spectr()")

#' @export
vis.immunr_spectr_nogene <- function(.data, ...) stop_v0_api("vis.immunr_spectr_nogene()")

#' @export
vis.immunr_top_prop <- function(.data, ...) stop_v0_api("vis.immunr_top_prop()")

#' @export
vis.immunr_tsne <- function(.data, ...) stop_v0_api("vis.immunr_tsne()")

#' @export
vis.step_failure_ignored <- function(.data, ...) stop_v0_api("vis.step_failure_ignored()")

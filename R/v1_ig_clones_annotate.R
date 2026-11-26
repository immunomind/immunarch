#' @title Analyse and annotate immunoglobulin receptors
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' A family of functions for analysing and adding biologically useful immunoglobulin
#' annotations to receptors. These annotations help connect receptor sequence
#' variation with clonal expansion and B-cell evolution.
#'
#' Each function has its own dedicated help page, so click on the function
#' names to open corresponding pages.
#'
#' ## Available functions
#'
#' 1) **Clonal family inference ([annotate_ig_clones()]).** Usually the very first step in
#'    any of the Ig analysis workflows. The function clusters Ig receptors to clonal families using
#'    heavy-chain or paired heavy-light junction similarity by Hamming distance.
#'
#' 2) **Germline annotation ([annotate_ig_germlines()]).** The first step in the
#'    somatic hypermutation (SHM) analysis workflows. This function annotates the
#'    Ig receptors with information about germline genes for the subsequent SHM analysis.
#'
#' 3) **Somatic hypermutation computation ([annotate_ig_shm()]).** The second and the final step in the
#'    SHM analysis workflows. This function computes the SHM per-receptor and aggregated
#'    statistics and annotates the receptors with them. After that you can use
#'    `airr_ig_shm()` to explore and analyse the repertoire-wide and family-wide statistics
#'    in a more convenient form.
#'
#' @seealso [annotate_ig_clones()], [annotate_ig_germlines()], [annotate_ig_shm()]
#' @name annotate_ig
#' @concept Immunoglobulin analysis
NULL


#' @title Build and annotate immunoglobulin clonal families
#'
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' Build clonal families from B-cell receptors and assign them to receptors.
#'
#' Receptors are compared within a subject or repertoire. A heavy-chain edge
#' requires a shared V/J union group, equal junction length, and normalized junction
#' Hamming distance at or below `heavy_threshold`. Paired mode independently
#' builds the corresponding light-chain graph and keeps receptors together only
#' when they share both heavy- and light-chain families.
#'
#' @param idata [immundata::ImmunData]. Paired or heavy-chain IG receptor data.
#' @param heavy_threshold,light_threshold Numeric scalar or data frame.
#'   Normalized Hamming-distance threshold strictly between `0` and `1`, or a
#'   complete table containing the `by_scope` column and a numeric `threshold`
#'   column.
#'   `light_threshold` defaults to `heavy_threshold`.
#' @param mode Character scalar. Annotation strategy: `"heavy"` uses only the
#'   heavy chain, `"paired"` combines independently inferred heavy- and
#'   light-chain families, and `"auto"` selects paired mode for a paired IG
#'   receptor schema and heavy mode otherwise.
#' @param by_scope Character scalar. Column that prevents clonal families from
#'   crossing biological subjects or repertoires. `"<auto>"` uses the first
#'   available subject-like column and otherwise the repertoire identifier.
#' @param junction_col,v_col,j_col Character scalar. Columns containing the
#'   nucleotide junction, V call, and J call. V and J calls are grouped at gene
#'   level, allowing overlapping ambiguous assignments.
#' @param heavy_locus Character scalar. Value identifying the heavy-chain
#'   locus. The default is `"IGH"`; all other loci declared by the receptor
#'   schema are treated as light chains in paired mode.
#' @param clone_col Character scalar. Name of the combined clone annotation.
#'   The default is `"clone_id"`.
#'
#' @details
#' Connected components use single linkage: when A connects to B and B connects
#' to C, all three receptors belong to one chain family even if A and C are
#' farther apart than the threshold. Hamming distance compares only junctions
#' of equal length, so it does not represent insertions or deletions.
#'
#' Junction nucleotides are converted to uppercase. Junctions containing
#' characters other than A, C, G, and T cannot create edges and remain isolated.
#'
#' For grouping, V/J calls are split on commas or `|`, whitespace and allele
#' suffixes are removed, and duplicate calls are ignored. Calls sharing both
#' a V candidate and a J candidate are connected transitively within each scope
#' and locus, before splitting by junction length or applying distance thresholds.
#' Missing or empty gene calls cannot connect groups. Original calls and receptor
#' IDs are preserved. The second example separately shows preprocessing when
#' normalized calls should also define receptor identity.
#'
#' @return [immundata::ImmunData]. The original rows and columns plus
#'   `clone_id` (or the name supplied through `clone_col`),
#'   `heavy_family_id`, and `light_family_id`. Both chain rows of a paired
#'   receptor receive the same labels. In heavy mode, `clone_id` equals
#'   `heavy_family_id` and `light_family_id` is missing.
#'
#' @examples
#' \dontshow{
#' # Release temporary DuckDB views from earlier examples after R resets the seed.
#' invisible(gc())
#' }
#' # Limit the number of threads used by the underlying DB for this example.
#' # Generally, you should NOT do this in your session.
#' db_exec("SET threads TO 1")
#'
#' # Six chain observations representing three paired B-cell receptors.
#' ig_chains <- tibble::tribble(
#'   ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
#'   "S1", "c1", "P1", "IGH", "AAAA", "IGHV1", "IGHJ1", 10L,
#'   "S1", "c1", "P1", "IGK", "CCCC", "IGKV1", "IGKJ1", 8L,
#'   "S1", "c2", "P1", "IGH", "AAAT", "IGHV1", "IGHJ1", 10L,
#'   "S1", "c2", "P1", "IGK", "CCCT", "IGKV1", "IGKJ1", 8L,
#'   "S1", "c3", "P1", "IGH", "GGGG", "IGHV2", "IGHJ2", 10L,
#'   "S1", "c3", "P1", "IGL", "TTTT", "IGLV1", "IGLJ1", 8L
#' )
#' input_file <- tempfile(fileext = ".tsv")
#' output_folder <- tempfile("ig-clones-")
#' readr::write_tsv(ig_chains, input_file)
#'
#' ig_schema <- make_receptor_schema(
#'   features = c("junction", "v_call", "j_call"),
#'   chains = c("IGH", "IGK|IGL")
#' )
#'
#' idata <- read_repertoires(
#'   path = input_file,
#'   schema = ig_schema,
#'   barcode_col = "cell",
#'   locus_col = "locus",
#'   umi_col = "umi",
#'   repertoire_schema = "sample",
#'   rename_columns = NULL,
#'   output_folder = output_folder,
#'   verbose = FALSE
#' )
#'
#' # c1 and c2 form one family; c3 remains separate.
#' annotated <- annotate_ig_clones(idata, heavy_threshold = 0.25)
#'
#' annotated
#'
#' \dontrun{
#' # Use heavy-chain similarity alone when the scientific question requires it.
#' annotate_ig_clones(
#'   idata,
#'   heavy_threshold = 0.25,
#'   mode = "heavy"
#' )
#'
#' # Canonicalize ambiguous V/J calls lazily before receptors are constructed.
#' # For example, "IGHV2*01, IGHV1*02" becomes "IGHV1|IGHV2".
#' canonicalize_vj <- function(data, ...) {
#'   data |>
#'     dplyr::mutate(dplyr::across(
#'       c(v_call, j_call),
#'       ~ dd$array_to_string(
#'         dd$list_sort(dd$list_distinct(dd$string_split(
#'           dd$regexp_replace(
#'             dd$regexp_replace(.x, "\\*[0-9]+", "", "g"),
#'             "\\s+", "", "g"
#'           ),
#'           ","
#'         ))),
#'         "|"
#'       )
#'     ))
#' }
#'
#' messy_chains <- ig_chains |>
#'   dplyr::mutate(
#'     v_call = dplyr::case_when(
#'       cell == "c2" & locus == "IGH" ~ " IGHV1*02, IGHV1*01 ",
#'       TRUE ~ paste0(v_call, "*01")
#'     ),
#'     j_call = paste0(j_call, "*01")
#'   )
#' messy_file <- tempfile(fileext = ".tsv")
#' canonical_folder <- tempfile("ig-clones-canonical-")
#' readr::write_tsv(messy_chains, messy_file)
#'
#' canonical_idata <- read_repertoires(
#'   path = messy_file,
#'   schema = ig_schema,
#'   barcode_col = "cell",
#'   locus_col = "locus",
#'   umi_col = "umi",
#'   repertoire_schema = "sample",
#'   rename_columns = NULL,
#'   preprocess = list(canonicalize_vj = canonicalize_vj),
#'   output_folder = canonical_folder,
#'   verbose = FALSE
#' )
#'
#' annotate_ig_clones(canonical_idata, heavy_threshold = 0.25)
#'
#' unlink(messy_file)
#' unlink(canonical_folder, recursive = TRUE)
#' }
#'
#' unlink(input_file)
#' unlink(output_folder, recursive = TRUE)
#' @name annotate_ig_clones
NULL

#' @keywords internal
annotate_ig_clones_impl <- function(
  idata,
  heavy_threshold,
  light_threshold = NULL,
  mode = c("auto", "heavy", "paired"),
  by_scope = "<auto>",
  junction_col = "junction",
  v_col = "v_call",
  j_col = "j_call",
  heavy_locus = "IGH",
  clone_col = "clone_id"
) {
  check_normalized_threshold <- function(threshold, argument) {
    if (is.null(threshold)) {
      cli::cli_abort("{.arg {argument}} must be a numeric scalar or a threshold table.")
    }

    if (is.data.frame(threshold) || !is.numeric(threshold) ||
      length(threshold) != 1L) {
      return(invisible())
    }

    if (is.na(threshold) || !is.finite(threshold)) {
      cli::cli_abort("{.arg {argument}} must be finite.")
    }

    if (threshold <= 0 || threshold >= 1) {
      cli::cli_abort("{.arg {argument}} must be strictly between 0 and 1.")
    }
  }

  mode <- match.arg(mode)
  checkmate::assert_string(by_scope)
  checkmate::assert_string(heavy_locus)
  checkmate::assert_string(clone_col)

  fixed_family_cols <- c("heavy_family_id", "light_family_id")
  if (clone_col %in% fixed_family_cols) {
    cli::cli_abort(
      "{.arg clone_col} must differ from {.field heavy_family_id} and {.field light_family_id}."
    )
  }

  annotation_cols <- colnames(idata$annotations)
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")
  locus_col <- immundata::imd_schema("locus")
  chains <- idata$schema_receptor$chains
  schema_loci <- unique(unlist(strsplit(chains, "|", fixed = TRUE)))
  light_loci <- setdiff(schema_loci, heavy_locus)
  paired_schema <- heavy_locus %in% schema_loci && length(light_loci) > 0L

  if (mode == "paired" && !paired_schema) {
    cli::cli_abort(
      "Paired mode requires a paired IGH and IGK/IGL receptor schema."
    )
  }

  if (mode == "auto") {
    mode <- if (paired_schema) "paired" else "heavy"
  }

  if (identical(by_scope, "<auto>")) {
    candidates <- intersect(
      c("subject_id", "donor", "donor_id", "patient_id", "patient", "subject"),
      annotation_cols
    )
    by_scope <- if (length(candidates)) {
      candidates[[1]]
    } else {
      repertoire_col
    }
  }

  if (!by_scope %in% annotation_cols) {
    cli::cli_abort(
      "Column {.field {by_scope}} selected by {.arg by_scope} is missing from {.field idata$annotations}."
    )
  }

  if (mode == "paired" && is.null(light_threshold)) {
    light_threshold <- heavy_threshold
  }

  check_normalized_threshold(heavy_threshold, "heavy_threshold")
  if (mode == "paired") {
    check_normalized_threshold(light_threshold, "light_threshold")
  }

  distance_idata <- idata |>
    dplyr::mutate(
      .ig_junction = dd$upper(!!rlang::sym(junction_col))
    ) |>
    dplyr::filter(dd$regexp_matches(
      !!rlang::sym(".ig_junction"),
      "^[ACGT]+$"
    ))
  distance_idata$schema_receptor$features <- unique(c(
    distance_idata$schema_receptor$features,
    ".ig_junction"
  ))

  families <- ig_clone_chain_components(
    idata,
    distance_idata,
    heavy_locus,
    heavy_threshold,
    "heavy_family_id",
    by_scope,
    locus_col,
    receptor_col,
    v_col,
    j_col
  )

  if (mode == "paired") {
    light <- ig_clone_chain_components(
      idata,
      distance_idata,
      light_loci,
      light_threshold,
      "light_family_id",
      by_scope,
      locus_col,
      receptor_col,
      v_col,
      j_col
    )

    families <- families |>
      dplyr::inner_join(light, by = c(by_scope, receptor_col))

    final_labels <- families |>
      dplyr::select(dplyr::all_of(c(
        by_scope,
        "heavy_family_id",
        "light_family_id"
      ))) |>
      dplyr::distinct() |>
      dplyr::arrange(!!!rlang::syms(c(
        by_scope,
        "heavy_family_id",
        "light_family_id"
      ))) |>
      dplyr::mutate(!!clone_col := dplyr::row_number())

    families <- families |>
      dplyr::left_join(
        final_labels,
        by = c(by_scope, "heavy_family_id", "light_family_id")
      )
  } else {
    families <- families |>
      dplyr::mutate(
        !!clone_col := !!rlang::sym("heavy_family_id"),
        light_family_id = NA_integer_
      )
  }

  join_by <- c(by_scope, receptor_col)
  names(join_by) <- join_by
  immundata::annotate(
    idata,
    annotations = families |>
      dplyr::select(dplyr::all_of(c(
        join_by,
        clone_col,
        fixed_family_cols
      ))),
    by = join_by,
    conflicts = "replace"
  )
}


ig_clone_chain_components <- function(
  idata,
  distance_idata,
  loci,
  threshold,
  output_col,
  by_scope,
  locus_col,
  receptor_col,
  v_col,
  j_col
) {
  chain_idata <- distance_idata |>
    dplyr::filter(!!rlang::sym(locus_col) %in% loci)
  chain_idata <- ig_clone_gene_groups(
    chain_idata,
    by = unique(c(locus_col, by_scope)),
    v_col = v_col,
    j_col = j_col
  )

  chain_edges <- dist_hamm(
    chain_idata,
    seq_col = ".ig_junction",
    by = unique(c(locus_col, by_scope, ".ig_vj_group")),
    max_dist = threshold,
    mode = "graph"
  )
  components <- clust_cc(
    idata,
    chain_edges,
    by = by_scope,
    cluster_col = ".ig_component",
    size_col = ".ig_component_size"
  )
  labels <- components |>
    dplyr::select(dplyr::all_of(c(by_scope, ".ig_component"))) |>
    dplyr::distinct() |>
    dplyr::arrange(!!!rlang::syms(c(by_scope, ".ig_component"))) |>
    dplyr::mutate(!!output_col := dplyr::row_number())

  components |>
    dplyr::left_join(labels, by = c(by_scope, ".ig_component")) |>
    dplyr::select(dplyr::all_of(c(by_scope, receptor_col, output_col)))
}


# Build the ambiguity graph on distinct call sets, not individual receptors.
ig_clone_gene_groups <- function(idata, by, v_col, j_col) {
  gene_calls <- function(column) {
    rlang::expr(dd$array_to_string(
      dd$list_sort(dd$list_distinct(dd$string_split(
        dd$regexp_replace(
          dd$replace(!!rlang::sym(column), "|", ","),
          "[[:space:]]|[*][^,]*", "", "g"
        ),
        ","
      ))),
      ","
    ))
  }

  grouped <- idata |>
    dplyr::mutate(
      .ig_v_genes = !!gene_calls(v_col),
      .ig_j_genes = !!gene_calls(j_col)
    ) |>
    dplyr::filter(
      dd$regexp_matches(.ig_v_genes, "[^,]"),
      dd$regexp_matches(.ig_j_genes, "[^,]")
    ) |>
    dplyr::mutate(
      # A length prefix makes the V/J key unambiguous without hashing.
      .ig_vj_key = dd$concat(
        dd$length(.ig_v_genes), ":", .ig_v_genes, .ig_j_genes
      )
    )
  source <- grouped$annotations |>
    dplyr::select(dplyr::all_of(c(
      by, ".ig_vj_key", ".ig_v_genes", ".ig_j_genes"
    ))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl()
  con <- dbplyr::remote_con(source)
  groups <- glue::glue_sql("{`by`*}", .con = con)

  # Each concrete V/J pair needs only a star of edges to connect its call sets.
  # Junction length deliberately plays no part in this preliminary graph.
  sql <- glue::glue_sql(
    "
    WITH calls AS ({source}), expanded AS (
      SELECT DISTINCT {groups}, \".ig_vj_key\", v.gene AS v, j.gene AS j
      FROM calls
      CROSS JOIN UNNEST(string_split(\".ig_v_genes\", ',')) AS v(gene)
      CROSS JOIN UNNEST(string_split(\".ig_j_genes\", ',')) AS j(gene)
      WHERE v.gene <> '' AND j.gene <> ''
    ), anchors AS (
      SELECT {groups}, \".ig_vj_key\" AS source,
        min(\".ig_vj_key\") OVER (PARTITION BY {groups}, v, j) AS target
      FROM expanded
    )
    SELECT DISTINCT {groups}, source, target
    FROM anchors
    WHERE source <> target",
    source = DBI::SQL(dbplyr::sql_render(source)),
    .con = con
  )

  edges <- dplyr::tbl(con, dbplyr::sql(sql)) |>
    duckplyr::as_duckdb_tibble()
  nodes <- grouped$annotations |>
    dplyr::select(dplyr::all_of(c(by, ".ig_vj_key"))) |>
    dplyr::distinct()
  labels <- clust_cc_tables(
    nodes, edges,
    from = "source", to = "target", by = by,
    node_col = ".ig_vj_key", cluster_col = ".ig_vj_group"
  )
  keys <- c(by, ".ig_vj_key")
  immundata::annotate(grouped, labels, by = stats::setNames(keys, keys))
}


#' @rdname annotate_ig_clones
#' @export
annotate_ig_clones <- register_method_annotate(
  core = annotate_ig_clones_impl,
  family = "annotate_ig",
  name = "clones",
  required_cols = c("junction_col", "v_col", "j_col"),
  need_repertoires = FALSE
)

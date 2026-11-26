make_dist_hamm_idata <- function() {
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")

  annotations <- tibble::tibble(
    !!receptor_col := c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
    !!repertoire_col := c(1L, 1L, 1L, 1L, 1L, 1L, 2L),
    cdr3_aa = c("AAA", "AAT", "ATT", "CCCC", "AAA", "AAA", "AAA"),
    v_call = c("V1", "V1", "V1", "V1", "V2", "V1", "V1"),
    j_call = c("J1", "J1", "J1", "J1", "J1", "J1", "J1")
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("cdr3_aa", "v_call", "j_call")
    ),
    annotations = annotations
  )
}


make_repeated_dist_hamm_idata <- function(n_identical) {
  receptor_col <- immundata::imd_schema("receptor")

  annotations <- tibble::tibble(
    !!receptor_col := seq_len(n_identical + 1L),
    cdr3_aa = c(rep("AAA", n_identical), "AAT"),
    v_call = "V1",
    j_call = "J1"
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("cdr3_aa", "v_call", "j_call")
    ),
    annotations = annotations
  )
}


make_paired_dist_hamm_idata <- function() {
  input_data <- tibble::tribble(
    ~barcode, ~locus, ~cdr3_aa, ~v_call, ~umis,
    "cell1", "TRA", "CAVR",  "TRAV1", 10L,
    "cell1", "TRB", "CASSA", "TRBV1", 20L,
    "cell2", "TRA", "CAVG",  "TRAV2", 15L,
    "cell2", "TRB", "CASSB", "TRBV2", 25L
  )

  schema <- immundata::make_receptor_schema(
    features = c("cdr3_aa", "v_call"),
    chains = c("TRA", "TRB")
  )

  path <- tempfile(fileext = ".tsv")
  readr::write_tsv(input_data, path)

  immundata::read_repertoires(
    path = path,
    schema = schema,
    barcode_col = "barcode",
    locus_col = "locus",
    umi_col = "umis",
    repertoire_schema = NULL,
    output_folder = tempfile("immundata-test-"),
    verbose = FALSE
  )
}


test_that("dist_hamm graph mode contracts identical sequences", {
  idata <- make_dist_hamm_idata()
  receptor_col <- immundata::imd_schema("receptor")

  graph <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1,
    mode = "graph"
  )

  expect_s3_class(graph, "immunarch_sequence_graph")
  expect_named(graph, c("edges", "membership"), ignore.order = FALSE)
  expect_s3_class(graph$edges, "duckplyr_df")
  expect_s3_class(graph$membership, "duckplyr_df")
  expect_named(
    graph$edges,
    c(
      "v_call", "j_call", "from", "to", "seq_len", "dist",
      "norm_dist", "sim"
    ),
    ignore.order = FALSE
  )
  expect_named(
    graph$membership,
    c("v_call", "j_call", receptor_col, "node_id"),
    ignore.order = FALSE
  )

  membership <- dplyr::collect(graph$membership)
  node_for <- function(receptor) {
    unique(membership$node_id[membership[[receptor_col]] == receptor])
  }

  expect_equal(nrow(membership), 7L)
  expect_equal(dplyr::n_distinct(membership$node_id), 5L)
  expect_length(node_for(1L), 1L)
  expect_equal(node_for(1L), node_for(6L))
  expect_equal(node_for(1L), node_for(7L))
  expect_false(identical(node_for(1L), node_for(2L)))
  expect_false(identical(node_for(1L), node_for(5L)))
})


test_that("dist_hamm graph mode returns unique-sequence distances", {
  idata <- make_dist_hamm_idata()
  receptor_col <- immundata::imd_schema("receptor")
  graph <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1,
    mode = "graph"
  )
  membership <- dplyr::collect(graph$membership)
  edges <- graph$edges |>
    dplyr::arrange(from, to) |>
    dplyr::collect()
  node_for <- function(receptor) {
    unique(membership$node_id[membership[[receptor_col]] == receptor])
  }
  edge_key <- function(from, to) {
    paste(sort(as.character(c(from, to))), collapse = ":")
  }

  actual_edges <- vapply(
    seq_len(nrow(edges)),
    function(i) edge_key(edges$from[[i]], edges$to[[i]]),
    character(1)
  )
  expected_edges <- c(
    edge_key(node_for(1L), node_for(2L)),
    edge_key(node_for(2L), node_for(3L))
  )

  expect_setequal(actual_edges, expected_edges)
  expect_true(all(edges$v_call == "V1"))
  expect_true(all(edges$j_call == "J1"))
  expect_true(all(edges$seq_len == 3L))
  expect_true(all(edges$dist == 1L))
  expect_true(all(edges$norm_dist == 1 / 3))
  expect_true(all(edges$sim == 1 - edges$norm_dist))

  incident_nodes <- unique(c(edges$from, edges$to))
  expect_false(node_for(4L) %in% incident_nodes)
  expect_false(node_for(5L) %in% incident_nodes)
})


test_that("dist_hamm graph contraction respects effective groups", {
  idata <- make_dist_hamm_idata()
  receptor_col <- immundata::imd_schema("receptor")
  repertoire_col <- immundata::imd_schema("repertoire")

  membership <- dist_hamm(
    idata,
    by = c(repertoire_col, "v_call", "j_call"),
    max_dist = 1,
    mode = "graph"
  )$membership |>
    dplyr::collect()
  node_for <- function(receptor) {
    unique(membership$node_id[membership[[receptor_col]] == receptor])
  }

  expect_equal(node_for(1L), node_for(6L))
  expect_false(identical(node_for(1L), node_for(7L)))
})


test_that("dist_hamm graph mode rejects reserved grouping columns", {
  receptor_col <- immundata::imd_schema("receptor")
  annotations <- tibble::tibble(
    !!receptor_col := 1:2,
    cdr3_aa = c("AAA", "AAT"),
    node_id = 1L,
    from = 1L,
    .idist_record_type = "input"
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "cdr3_aa"),
    annotations = annotations
  )

  for (group_col in c("node_id", "from", ".idist_record_type")) {
    expect_error(
      dist_hamm(idata, by = group_col, mode = "graph"),
      "reserved.*output"
    )
  }
})


test_that("dist_hamm preserves duplicate sequence-derived receptor pairs", {
  receptor_col <- immundata::imd_schema("receptor")
  annotations <- tibble::tibble(
    !!receptor_col := c(1L, 1L, 2L, 2L),
    cdr3_aa = c("AAA", "AAT", "AAC", "AGT")
  ) |>
    duckplyr::as_duckdb_tibble(prudence = "stingy")
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(features = "cdr3_aa"),
    annotations = annotations
  )

  out <- dist_hamm(idata) |>
    dplyr::arrange(dist) |>
    dplyr::collect()

  expect_equal(out$imd_receptor_id_1, rep(1L, 4L))
  expect_equal(out$imd_receptor_id_2, rep(2L, 4L))
  expect_equal(out$dist, c(1, 1, 1, 2))
})


test_that("dist_hamm graph size is independent of duplicate multiplicity", {
  graph_sizes <- lapply(c(2L, 20L), function(n_identical) {
    graph <- dist_hamm(
      make_repeated_dist_hamm_idata(n_identical),
      by = c("v_call", "j_call"),
      max_dist = 1,
      mode = "graph"
    )

    c(
      nodes = graph$membership |>
        dplyr::distinct(node_id) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::pull(n),
      edges = graph$edges |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::pull(n),
      memberships = graph$membership |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::pull(n)
    )
  })

  expect_equal(graph_sizes[[1]][["nodes"]], 2)
  expect_equal(graph_sizes[[2]][["nodes"]], 2)
  expect_equal(graph_sizes[[1]][["edges"]], 1)
  expect_equal(graph_sizes[[2]][["edges"]], 1)
  expect_equal(graph_sizes[[1]][["memberships"]], 3)
  expect_equal(graph_sizes[[2]][["memberships"]], 21)
})


test_that("dist_hamm returns lazy receptor-level upper-triangle edges", {
  idata <- make_dist_hamm_idata()

  lazy <- dist_hamm(
    idata,
    by = c("v_call", "j_call")
  )

  expect_s3_class(lazy, "duckplyr_df")

  invisible(gc())

  out <- lazy |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  expect_named(
    out,
    c(
      "v_call", "j_call", "imd_receptor_id_1", "imd_receptor_id_2",
      "seq_len", "dist", "norm_dist", "sim"
    ),
    ignore.order = FALSE
  )
  expect_equal(
    out$imd_receptor_id_1,
    c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 6L)
  )
  expect_equal(
    out$imd_receptor_id_2,
    c(2L, 3L, 6L, 7L, 3L, 6L, 7L, 6L, 7L, 7L)
  )
  expect_equal(out$dist, c(1, 2, 0, 0, 1, 1, 1, 2, 2, 0))
  expect_equal(out$norm_dist, out$dist / 3)
  expect_equal(out$sim, 1 - out$norm_dist)
  expect_true(all(out$seq_len == 3L))
})


test_that("dist_hamm materialization survives an RNG reset", {
  idata <- make_dist_hamm_idata()

  withr::local_preserve_seed()
  set.seed(431)
  first <- dist_hamm(idata)
  invisible(gc())

  set.seed(431)
  second <- dist_hamm(idata)

  arrange_edges <- function(edges) {
    edges |>
      dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
      dplyr::collect()
  }
  expect_equal(arrange_edges(second), arrange_edges(first))
})


test_that("dist_hamm computes paired-chain distances separately by locus", {
  idata <- make_paired_dist_hamm_idata()

  out <- dist_hamm(idata) |>
    dplyr::arrange(locus) |>
    dplyr::collect()

  expect_named(
    out,
    c(
      "locus", "imd_receptor_id_1", "imd_receptor_id_2",
      "seq_len", "dist", "norm_dist", "sim"
    ),
    ignore.order = FALSE
  )
  expect_equal(out$locus, c("TRA", "TRB"))
  expect_equal(out$imd_receptor_id_1, c(1L, 1L))
  expect_equal(out$imd_receptor_id_2, c(2L, 2L))
  expect_equal(out$seq_len, c(4L, 5L))
  expect_equal(out$dist, c(1, 1))
  expect_equal(out$norm_dist, c(0.25, 0.2))

  bounded <- dist_hamm(
    idata,
    by = "locus",
    max_dist = 0.21
  ) |>
    dplyr::collect()

  expect_equal(bounded$locus, "TRB")
  expect_equal(bounded$dist, 1)
  expect_equal(bounded$norm_dist, 0.2)
})


test_that("dist_hamm and locus filtering are order independent", {
  idata <- make_paired_dist_hamm_idata()

  distance_then_filter <- dist_hamm(idata) |>
    dplyr::filter(.data$locus == "TRA") |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  filter_then_distance <- idata |>
    dplyr::filter(.data$locus == "TRA") |>
    dist_hamm() |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(distance_then_filter, filter_then_distance)
})


test_that("dist_hamm always returns normalized Hamming distance", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = c("v_call", "j_call")
  ) |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(
    out$norm_dist,
    c(1 / 3, 2 / 3, 0, 0, 1 / 3, 1 / 3, 1 / 3, 2 / 3, 2 / 3, 0)
  )
})


test_that("dist_hamm bounded modes use exact candidate generation", {
  idata <- make_dist_hamm_idata()

  raw <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1
  ) |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  normalized <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1 / 3
  ) |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  similar <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    min_sim = 2 / 3
  ) |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(
    raw$imd_receptor_id_1,
    c(1L, 1L, 1L, 2L, 2L, 2L, 6L)
  )
  expect_equal(
    raw$imd_receptor_id_2,
    c(2L, 6L, 7L, 3L, 6L, 7L, 7L)
  )
  expect_equal(raw$dist, c(1, 0, 0, 1, 1, 1, 0))
  expect_equal(normalized$imd_receptor_id_1, raw$imd_receptor_id_1)
  expect_equal(normalized$imd_receptor_id_2, raw$imd_receptor_id_2)
  expect_equal(normalized$dist, raw$dist)
  expect_equal(normalized$norm_dist, raw$dist / 3)
  expect_equal(similar, normalized)
})


test_that("dist_hamm applies keyed normalized distance thresholds", {
  thresholds <- tibble::tribble(
    ~subject_id, ~threshold,
    "P1", 0.25,
    "P2", 0.50
  )

  out <- dist_hamm(
    make_grouped_distance_test_idata(),
    by = c("subject_id", "v_call", "j_call"),
    max_dist = thresholds
  ) |>
    dplyr::collect()

  expect_equal(as.integer(table(out$subject_id)), c(3L, 5L))
  expect_equal(max(out$norm_dist[out$subject_id == "P1"]), 0.25)
  expect_equal(max(out$norm_dist[out$subject_id == "P2"]), 0.50)
  expect_false("threshold" %in% names(out))
})


test_that("dist_hamm supports multi-column threshold keys", {
  thresholds <- tibble::tribble(
    ~subject_id, ~v_call, ~threshold,
    "P1", "V1", 0.25,
    "P2", "V1", 0.50
  )

  out <- dist_hamm(
    make_grouped_distance_test_idata(),
    by = c("subject_id", "v_call", "j_call"),
    max_dist = thresholds
  ) |>
    dplyr::collect()

  expect_equal(as.integer(table(out$subject_id)), c(3L, 5L))
})


test_that("dist_hamm similarity bounds include identity", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(idata, min_sim = 1) |>
    dplyr::collect()

  expect_true(all(out$dist == 0))
  expect_true(all(out$sim == 1))
})


test_that("dist_hamm grouping is optional and length matching is automatic", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = NULL,
    max_dist = 1e-6
  ) |>
    dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
    dplyr::collect()

  expect_equal(out$imd_receptor_id_1, c(1L, 1L, 1L, 5L, 5L, 6L))
  expect_equal(out$imd_receptor_id_2, c(5L, 6L, 7L, 6L, 7L, 7L))
  expect_true(all(out$dist == 0))
  expect_false(any(
    out$imd_receptor_id_1 == 4L | out$imd_receptor_id_2 == 4L
  ))
})


test_that("dist_hamm uses annotation columns in by", {
  idata <- make_dist_hamm_idata()
  repertoire_col <- immundata::imd_schema("repertoire")

  global <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 1e-6
  ) |>
    dplyr::collect()

  within_repertoire <- dist_hamm(
    idata,
    by = c(repertoire_col, "v_call", "j_call"),
    max_dist = 1e-6
  ) |>
    dplyr::collect()

  expect_true(any(
    global$imd_receptor_id_1 == 1L & global$imd_receptor_id_2 == 7L
  ))
  expect_equal(within_repertoire$imd_receptor_id_1, 1L)
  expect_equal(within_repertoire$imd_receptor_id_2, 6L)
  expect_equal(within_repertoire[[repertoire_col]], 1L)
})


test_that("dist_hamm samples receptors within an inferred subject column", {
  out <- dist_hamm(
    make_grouped_distance_test_idata(),
    by = c("subject_id", "v_call", "j_call"),
    sample_n = 2L
  ) |>
    dplyr::collect()

  expect_equal(sort(out$subject_id), c("P1", "P2"))
  expect_equal(nrow(out), 2L)
  expect_equal(
    length(unique(c(
      out$imd_receptor_id_1,
      out$imd_receptor_id_2
    ))),
    4L
  )
})


test_that("dist_hamm supports global receptor sampling", {
  out <- dist_hamm(
    make_grouped_distance_test_idata(),
    by = c("v_call", "j_call"),
    sample_n = 3L,
    sample_by = NULL
  ) |>
    dplyr::collect()

  expect_equal(nrow(out), 3L)
  expect_equal(
    length(unique(c(
      out$imd_receptor_id_1,
      out$imd_receptor_id_2
    ))),
    3L
  )
})


test_that("dist_hamm sampling also works with bounded distances", {
  out <- dist_hamm(
    make_grouped_distance_test_idata(),
    by = c("subject_id", "v_call", "j_call"),
    max_dist = 4L,
    sample_n = 2L
  ) |>
    dplyr::collect()

  expect_equal(sort(out$subject_id), c("P1", "P2"))
  expect_equal(nrow(out), 2L)
})


test_that("dist_hamm block search exactly matches filtered direct distance", {
  idata <- make_dist_hamm_idata()

  direct_raw <- dist_hamm(
    idata,
    by = c("v_call", "j_call")
  ) |>
    dplyr::collect()

  for (bound in c(1, 2, 10)) {
    expected <- direct_raw |>
      dplyr::filter(.data$dist <= .env$bound) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2)

    actual <- dist_hamm(
      idata,
      by = c("v_call", "j_call"),
      max_dist = bound
    ) |>
      dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
      dplyr::collect()

    expect_equal(actual, expected)
  }

  for (bound in c(1e-6, 1 / 3, 2 / 3)) {
    expected <- direct_raw |>
      dplyr::filter(.data$norm_dist <= .env$bound + 1e-12) |>
      dplyr::arrange(.data$imd_receptor_id_1, .data$imd_receptor_id_2)

    actual <- dist_hamm(
      idata,
      by = c("v_call", "j_call"),
      max_dist = bound
    ) |>
      dplyr::arrange(imd_receptor_id_1, imd_receptor_id_2) |>
      dplyr::collect()

    expect_equal(actual, expected, tolerance = 1e-12)
  }
})


test_that("dist_hamm does not require repertoire aggregation", {
  idata <- make_dist_hamm_idata()

  out <- dist_hamm(
    idata,
    by = c("v_call", "j_call"),
    max_dist = 0
  )

  expect_s3_class(out, "duckplyr_df")
  expect_null(idata$repertoires)
  expect_false(immundata::imd_schema("repertoire") %in% colnames(out))
})


test_that("dist_hamm validates columns and bounds", {
  idata <- make_dist_hamm_idata()

  expect_error(dist_hamm(idata, mode = "unknown"), "mode")
  expect_error(dist_hamm(idata, seq_col = "missing"), "missing")
  expect_error(dist_hamm(idata, by = "missing"), "missing")
  expect_error(dist_hamm(idata, max_dist = 1.5), "integer")
  expect_error(dist_hamm(idata, min_sim = 1.1), "<= 1")
  expect_error(
    dist_hamm(idata, max_dist = 0.1, min_sim = 0.9),
    "only one"
  )
  expect_error(dist_hamm(idata, sample_n = -1L), ">= 0", fixed = TRUE)
  expect_error(
    dist_hamm(idata, sample_n = 2L),
    "absent from.*by"
  )
  expect_error(
    dist_hamm(
      idata,
      by = c("v_call", "j_call"),
      sample_n = 2L,
      sample_by = "subject_id"
    ),
    "missing"
  )
  expect_no_error(
    dist_hamm(idata, sample_n = 0L, sample_by = "missing")
  )

  idata_without_sequence_schema <- immundata::ImmunData$new(
    schema = "v_call",
    annotations = idata$annotations
  )
  expect_error(
    dist_hamm(idata_without_sequence_schema, seq_col = "cdr3_aa"),
    "receptor schema"
  )
})


test_that("dist_hamm validates keyed thresholds", {
  idata <- make_grouped_distance_test_idata()
  by <- c("subject_id", "v_call", "j_call")
  valid <- tibble::tibble(
    subject_id = c("P1", "P2"),
    threshold = c(0.25, 0.50)
  )

  expect_error(
    dist_hamm(idata, by = by, max_dist = valid[1, ]),
    "cover every input"
  )
  expect_error(
    dist_hamm(
      idata,
      by = by,
      max_dist = dplyr::bind_rows(
        valid,
        tibble::tibble(subject_id = "P3", threshold = 0.25)
      )
    ),
    "absent from the input"
  )
  expect_error(
    dist_hamm(idata, by = by, max_dist = dplyr::bind_rows(valid, valid[1, ])),
    "must be unique"
  )
  expect_error(
    dist_hamm(idata, by = by, max_dist = dplyr::mutate(valid, threshold = 0)),
    "between 0 and 1"
  )
  expect_error(
    dist_hamm(idata, by = c("v_call", "j_call"), max_dist = valid),
    "included in.*by"
  )
  expect_error(
    dist_hamm(idata, by = by, max_dist = valid, min_sim = 0.5),
    "only one"
  )
})


test_that("dist_hamm is registered in the dist family", {
  expect_true("hamm" %in% ls(IMMUNARCH_METHOD_REGISTRY[["dist"]]))
})

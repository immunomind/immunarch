# Connected-component clustering tests for contracted sequence graphs.


make_clust_cc_idata <- function(rows, chains = "IGH") {
  path <- tempfile(fileext = ".tsv")
  readr::write_tsv(rows, path)

  immundata::read_repertoires(
    path = path,
    schema = immundata::make_receptor_schema(
      features = c("sequence", "comparison"),
      chains = chains
    ),
    barcode_col = "cell",
    locus_col = "locus",
    umi_col = "umi",
    repertoire_schema = "sample",
    rename_columns = NULL,
    output_folder = tempfile("immundata-test-"),
    verbose = FALSE
  )
}


clust_cc_rows <- function(
  cell,
  sequence,
  scope = "A",
  comparison = "same",
  sample = "S1",
  locus = "IGH"
) {
  tibble::tibble(
    sample = sample,
    cell = cell,
    scope = scope,
    comparison = comparison,
    locus = locus,
    sequence = sequence,
    umi = 1L
  )
}


clust_cc_graph <- function(idata, by = "comparison", max_dist = 1) {
  dist_hamm(
    idata,
    seq_col = "sequence",
    by = by,
    max_dist = max_dist,
    mode = "graph"
  )
}


collect_clust_cc <- function(result, by = "scope") {
  receptor_col <- immundata::imd_schema("receptor")

  result |>
    dplyr::arrange(!!!rlang::syms(c(by, receptor_col))) |>
    dplyr::collect() |>
    tibble::as_tibble()
}


clust_cc_receptors <- function(idata) {
  receptor_col <- immundata::imd_schema("receptor")

  idata$annotations |>
    dplyr::select(dplyr::all_of(c("scope", "cell", receptor_col))) |>
    dplyr::distinct() |>
    dplyr::collect()
}


expect_clust_cc_partition <- function(
  actual,
  receptors,
  expected_components,
  cluster_col = "cluster_id",
  size_col = "cluster_size"
) {
  receptor_col <- immundata::imd_schema("receptor")
  expected <- dplyr::bind_rows(lapply(
    expected_components,
    function(cells) {
      members <- receptors |>
        dplyr::filter(cell %in% cells) |>
        dplyr::distinct(scope, !!rlang::sym(receptor_col))
      members[[cluster_col]] <- min(members[[receptor_col]])
      members[[size_col]] <- dplyr::n_distinct(members[[receptor_col]])
      members
    }
  )) |>
    dplyr::select(dplyr::all_of(c(
      "scope", receptor_col, cluster_col, size_col
    ))) |>
    dplyr::arrange(scope, !!rlang::sym(receptor_col))

  actual <- actual |>
    dplyr::select(dplyr::all_of(names(expected))) |>
    dplyr::arrange(scope, !!rlang::sym(receptor_col))

  testthat::expect_equal(actual, expected)
}


testthat::test_that("clust_cc rejects ordinary receptor edge tables", {
  idata <- make_clust_cc_idata(clust_cc_rows(
    cell = c("a", "b"),
    sequence = c("AAAA", "AAAT")
  ))
  receptor_col <- immundata::imd_schema("receptor")
  ordinary_edges <- idata$annotations |>
    dplyr::select(dplyr::all_of(c("scope", receptor_col))) |>
    dplyr::distinct() |>
    duckplyr::as_tbl() |>
    dplyr::arrange(!!rlang::sym(receptor_col)) |>
    dplyr::mutate(
      imd_receptor_id_2 = dplyr::lead(!!rlang::sym(receptor_col))
    ) |>
    dplyr::rename(
      imd_receptor_id_1 = !!rlang::sym(receptor_col)
    ) |>
    dplyr::filter(!is.na(imd_receptor_id_2))

  testthat::expect_error(
    clust_cc(idata, ordinary_edges, by = "scope"),
    "mode = [\"']graph[\"']"
  )
})


testthat::test_that(
  "clust_cc labels receptor components and retains isolated receptors",
  {
    rows <- clust_cc_rows(
      cell = c("a", "b", "c", "duplicate", "isolated"),
      sequence = c("AAAA", "AAAT", "AATT", "AAAA", "GGGG")
    )
    idata <- make_clust_cc_idata(rows)
    receptors <- clust_cc_receptors(idata)

    result <- clust_cc(
      idata,
      clust_cc_graph(idata),
      by = "scope",
      cluster_col = "component",
      size_col = "receptor_count"
    )
    actual <- collect_clust_cc(result)

    testthat::expect_s3_class(result, "immunarch_res_clust_cc")
    testthat::expect_named(
      actual,
      c("scope", immundata::imd_schema("receptor"), "component", "receptor_count")
    )
    expect_clust_cc_partition(
      actual,
      receptors,
      list(
        c("a", "b", "c", "duplicate"),
        "isolated"
      ),
      cluster_col = "component",
      size_col = "receptor_count"
    )
  }
)


testthat::test_that("receptors absent from graph membership remain singletons", {
  receptor_col <- immundata::imd_schema("receptor")
  idata <- make_clust_cc_idata(clust_cc_rows(
    cell = c("a", "b", "omitted"),
    sequence = c("AAAA", "AAAT", "GGGG")
  ))
  graph <- idata |>
    dplyr::filter(cell != "omitted") |>
    clust_cc_graph(by = NULL)
  omitted_id <- idata$annotations |>
    dplyr::filter(cell == "omitted") |>
    dplyr::select(dplyr::all_of(receptor_col)) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::pull(dplyr::all_of(receptor_col))

  omitted <- clust_cc(idata, graph, by = NULL) |>
    dplyr::filter(!!rlang::sym(receptor_col) == omitted_id) |>
    dplyr::collect()

  testthat::expect_equal(omitted$cluster_id, omitted_id)
  testthat::expect_equal(omitted$cluster_size, 1)
})


testthat::test_that(
  "one receptor bridges multiple sequence signatures without receptor pairs",
  {
    rows <- dplyr::bind_rows(
      clust_cc_rows("left", "AAAA", locus = "IGH"),
      clust_cc_rows("left", "GGGG", locus = "IGK"),
      clust_cc_rows("bridge", "AAAT", locus = "IGH"),
      clust_cc_rows("bridge", "CCCC", locus = "IGK"),
      clust_cc_rows("right", "TTTT", locus = "IGH"),
      clust_cc_rows("right", "CCCT", locus = "IGK")
    )
    idata <- make_clust_cc_idata(rows, chains = c("IGH", "IGK"))
    receptors <- clust_cc_receptors(idata)
    graph <- clust_cc_graph(idata, by = "comparison")

    membership_counts <- graph$membership |>
      dplyr::summarise(n = dplyr::n(), .by = immundata::imd_schema("receptor")) |>
      dplyr::collect()
    testthat::expect_true(all(membership_counts$n == 2L))

    actual <- clust_cc(idata, graph, by = "scope") |>
      collect_clust_cc()

    testthat::expect_equal(nrow(actual), 3L)
    expect_clust_cc_partition(
      actual,
      receptors,
      list(c("left", "bridge", "right"))
    )
    testthat::expect_true(all(actual$cluster_size == 3L))
  }
)


testthat::test_that("graph membership remains scoped by locus", {
  rows <- dplyr::bind_rows(
    clust_cc_rows("left", "AAAA", locus = "IGH"),
    clust_cc_rows("left", "GGGG", locus = "IGK"),
    clust_cc_rows("bridge", "AAAT", locus = "IGH"),
    clust_cc_rows("bridge", "CCCC", locus = "IGK"),
    clust_cc_rows("right", "TTTT", locus = "IGH"),
    clust_cc_rows("right", "CCCT", locus = "IGK")
  )
  idata <- make_clust_cc_idata(rows, chains = c("IGH", "IGK"))
  receptor_col <- immundata::imd_schema("receptor")
  receptors <- idata$annotations |>
    dplyr::select("locus", "cell", dplyr::all_of(receptor_col)) |>
    dplyr::distinct() |>
    dplyr::collect()

  actual <- clust_cc(
    idata,
    clust_cc_graph(idata, by = "locus"),
    by = "locus"
  ) |>
    collect_clust_cc(by = "locus") |>
    dplyr::left_join(receptors, by = c("locus", receptor_col))

  cluster <- function(wanted_locus, wanted_cell) {
    actual |>
      dplyr::filter(
        locus == .env$wanted_locus,
        cell == .env$wanted_cell
      ) |>
      dplyr::pull("cluster_id")
  }

  testthat::expect_equal(cluster("IGH", "left"), cluster("IGH", "bridge"))
  testthat::expect_false(cluster("IGH", "left") == cluster("IGH", "right"))
  testthat::expect_equal(cluster("IGK", "bridge"), cluster("IGK", "right"))
  testthat::expect_false(cluster("IGK", "left") == cluster("IGK", "bridge"))
})


testthat::test_that(
  "comparison groups do not replace clustering scopes",
  {
    rows <- clust_cc_rows(
      cell = c("a1", "a2", "b1", "b2", "null1", "null2"),
      sequence = c("AAAA", "AAAT", "AAAA", "AATT", "CCCC", "CCCT"),
      scope = c("A", "A", "B", "B", NA, NA),
      comparison = "global"
    )
    idata <- make_clust_cc_idata(rows)
    receptors <- clust_cc_receptors(idata)
    graph <- clust_cc_graph(idata, by = "comparison")

    actual <- clust_cc(idata, graph, by = "scope") |>
      collect_clust_cc()

    expect_clust_cc_partition(
      actual,
      receptors,
      list(c("a1", "a2"), "b1", "b2", c("null1", "null2"))
    )
  }
)


testthat::test_that("independent null scopes remain independent", {
  rows <- clust_cc_rows(
    cell = c("a1", "a2", "null1", "null2"),
    sequence = c("AAAA", "AAAT", "AAAA", "AAAT"),
    scope = c("A", "A", NA, NA),
    comparison = "global"
  )
  idata <- make_clust_cc_idata(rows)
  receptors <- clust_cc_receptors(idata)

  actual <- clust_cc(idata, clust_cc_graph(idata), by = "scope") |>
    collect_clust_cc()

  expect_clust_cc_partition(
    actual,
    receptors,
    list(c("a1", "a2"), c("null1", "null2"))
  )
})


testthat::test_that("transitive paths, cycles, and graph row order are stable", {
  chain_rows <- clust_cc_rows(
    cell = c("chain1", "chain2", "chain3"),
    sequence = c("AAAA", "AAAT", "AATT"),
    comparison = "chain"
  )
  cycle_rows <- clust_cc_rows(
    cell = c("cycle1", "cycle2", "cycle3"),
    sequence = c("CCCC", "CCCA", "CCCT"),
    comparison = "cycle"
  )
  idata <- make_clust_cc_idata(dplyr::bind_rows(chain_rows, cycle_rows))
  graph <- clust_cc_graph(idata)
  shuffled <- list(
    edges = graph$edges |>
      dplyr::arrange(dplyr::desc(from), dplyr::desc(to)),
    membership = graph$membership |>
      dplyr::arrange(dplyr::desc(node_id))
  )
  class(shuffled) <- class(graph)

  original <- clust_cc(idata, graph, by = "scope") |>
    collect_clust_cc()
  reordered <- clust_cc(idata, shuffled, by = "scope") |>
    collect_clust_cc()

  testthat::expect_equal(reordered, original)
  expect_clust_cc_partition(
    original,
    clust_cc_receptors(idata),
    list(
      c("chain1", "chain2", "chain3"),
      c("cycle1", "cycle2", "cycle3")
    )
  )
})


testthat::test_that("empty receptor universes return an empty typed result", {
  idata <- make_clust_cc_idata(clust_cc_rows("seed", "AAAA")) |>
    dplyr::filter(FALSE)
  graph <- clust_cc_graph(idata)

  result <- clust_cc(
    idata,
    graph,
    by = "scope",
    cluster_col = "component",
    size_col = "receptor_count"
  )
  actual <- dplyr::collect(result)

  testthat::expect_s3_class(result, "immunarch_res_clust_cc")
  testthat::expect_named(
    actual,
    c("scope", immundata::imd_schema("receptor"), "component", "receptor_count")
  )
  testthat::expect_equal(nrow(actual), 0L)
})

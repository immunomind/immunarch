make_ig_clone_idata <- function(
  rows,
  chains = c("IGH", "IGK|IGL"),
  features = c("junction", "v_call", "j_call")
) {
  path <- tempfile(fileext = ".tsv")
  readr::write_tsv(rows, path)

  idata <- immundata::read_repertoires(
    path = path,
    schema = immundata::make_receptor_schema(
      features = features,
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


paired_ig_rows <- function(cells) {
  heavy <- cells |>
    dplyr::transmute(
      .data$sample,
      .data$cell,
      .data$subject,
      locus = "IGH",
      junction = .data$h_junction,
      v_call = .data$h_v_call,
      j_call = .data$h_j_call,
      umi = 10L
    )

  light <- cells |>
    dplyr::transmute(
      .data$sample,
      .data$cell,
      .data$subject,
      locus = .data$l_locus,
      junction = .data$l_junction,
      v_call = .data$l_v_call,
      j_call = .data$l_j_call,
      umi = 8L
    )

  dplyr::bind_rows(heavy, light)
}


family_partition <- function(data, member_col, family_col = "family_id") {
  collected <- data |>
    dplyr::select(dplyr::all_of(c(member_col, family_col))) |>
    dplyr::distinct() |>
    dplyr::collect()

  split(collected[[member_col]], collected[[family_col]]) |>
    lapply(sort) |>
    vapply(paste, collapse = ",", character(1)) |>
    unname() |>
    sort()
}


collect_cell_families <- function(idata) {
  idata$annotations |>
    dplyr::select(
      dplyr::all_of(c(
        "cell",
        "locus",
        "family_id",
        "heavy_family_id",
        "light_family_id"
      ))
    ) |>
    dplyr::collect()
}


testthat::test_that("paired families use heavy and light junction distances", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AATT", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTT", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1,c2", "c3"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    "c1,c2,c3"
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "light_family_id"),
    sort(c("c1,c2", "c3"))
  )
})


testthat::test_that("paired families intersect chain components, not edges", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AATT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCA", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    light_threshold = 0.25
  )

  actual <- collect_cell_families(out)

  # Heavy edges are c1-c2 and c2-c3. Light edges are c1-c2 and c1-c3.
  # Both chain graphs have one component, despite having different edge sets.
  testthat::expect_equal(family_partition(actual, "cell"), "c1,c2,c3")
  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    "c1,c2,c3"
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "light_family_id"),
    "c1,c2,c3"
  )
})


testthat::test_that("light threshold is independent of heavy threshold", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    light_threshold = 0.20
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    "c1,c2"
  )
  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1", "c2"))
  )
})


testthat::test_that("IGK and IGL define incompatible light families", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGL", "CCCC", "IGLV1*01", "IGLJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    "c1,c2"
  )
  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1", "c2"))
  )
})


testthat::test_that("different heavy or light junction lengths cannot join", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCCC", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )
  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    sort(c("c1,c3", "c2"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "light_family_id"),
    sort(c("c1,c2", "c3"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1", "c2", "c3"))
  )
})


testthat::test_that("disjoint heavy or light J genes cannot join", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ2*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ2*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )
  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    sort(c("c1,c3", "c2"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "light_family_id"),
    sort(c("c1,c2", "c3"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1", "c2", "c3"))
  )
})


testthat::test_that("lowercase junctions normalize and ambiguous bases stay singleton", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "aaaa", "IGHV1*01", "IGHJ1*01",
    "IGK", "cccc", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AAAN", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCN", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    light_threshold = 0.25
  )

  testthat::expect_equal(
    family_partition(out$annotations, "cell"),
    sort(c("c1,c2", "c3"))
  )
})


testthat::test_that("automatic scope joins samples within subjects only", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S2", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "B", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    by_scope = "<auto>"
  )

  testthat::expect_equal(
    family_partition(out$annotations, "cell"),
    sort(c("c1,c2", "c3"))
  )
})


testthat::test_that("automatic scope falls back to repertoires", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "unused", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S2", "c2", "unused", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )
  rows <- paired_ig_rows(cells) |>
    dplyr::select(-"subject")

  out <- annotate_ig_clones(
    make_ig_clone_idata(rows),
    heavy_threshold = 0.25,
    by_scope = "<auto>"
  )

  testthat::expect_equal(
    family_partition(out$annotations, "cell"),
    sort(c("c1", "c2"))
  )
})


testthat::test_that("family partitions are stable after shuffling input rows", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AATT", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTT", "IGKV1*01", "IGKJ1*01",
    "S1", "c4", "A", "GGGG", "IGHV3*01", "IGHJ3*01",
    "IGL", "GGGG", "IGLV2*01", "IGLJ2*01"
  )
  rows <- paired_ig_rows(cells)
  shuffled_rows <- rows[rev(seq_len(nrow(rows))), ]

  original <- annotate_ig_clones(
    make_ig_clone_idata(rows),
    heavy_threshold = 0.25
  )
  shuffled <- annotate_ig_clones(
    make_ig_clone_idata(shuffled_rows),
    heavy_threshold = 0.25
  )

  for (family_col in c(
    "family_id",
    "heavy_family_id",
    "light_family_id"
  )) {
    testthat::expect_equal(
      family_partition(original$annotations, "cell", family_col),
      family_partition(shuffled$annotations, "cell", family_col)
    )
  }
})


testthat::test_that("heavy and light thresholds accept scope-keyed tables", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "a1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "a2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "b1", "B", "GGGG", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTT", "IGKV1*01", "IGKJ1*01",
    "S1", "b2", "B", "GGGA", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTA", "IGKV1*01", "IGKJ1*01"
  )
  heavy_thresholds <- tibble::tribble(
    ~subject, ~threshold,
    "A", 0.25,
    "B", 0.20
  )
  light_thresholds <- tibble::tribble(
    ~subject, ~threshold,
    "A", 0.25,
    "B", 0.25
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = heavy_thresholds,
    light_threshold = light_thresholds,
    by_scope = "subject"
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("a1,a2", "b1", "b2"))
  )
  testthat::expect_equal(
    family_partition(actual, "cell", "light_family_id"),
    sort(c("a1,a2", "b1,b2"))
  )
})


testthat::test_that("scalar thresholds must be finite normalized distances", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )
  idata <- make_ig_clone_idata(paired_ig_rows(cells))
  out_of_range <- c(-0.1, 0, 1, 1.1)
  non_finite <- c(NA_real_, NaN, Inf)

  for (threshold in out_of_range) {
    testthat::expect_error(
      annotate_ig_clones(
        idata,
        heavy_threshold = threshold
      ),
      "heavy_threshold.*strictly between 0 and 1"
    )
    testthat::expect_error(
      annotate_ig_clones(
        idata,
        heavy_threshold = 0.25,
        light_threshold = threshold
      ),
      "light_threshold.*strictly between 0 and 1"
    )
  }
  for (threshold in non_finite) {
    testthat::expect_error(
      annotate_ig_clones(
        idata,
        heavy_threshold = threshold
      ),
      "heavy_threshold.*finite"
    )
    testthat::expect_error(
      annotate_ig_clones(
        idata,
        heavy_threshold = 0.25,
        light_threshold = threshold
      ),
      "light_threshold.*finite"
    )
  }
})


testthat::test_that("threshold tables require complete unique valid scopes", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "a1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "b1", "B", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01"
  )
  idata <- make_ig_clone_idata(paired_ig_rows(cells))

  testthat::expect_error(
    annotate_ig_clones(
      idata,
      heavy_threshold = tibble::tibble(
        subject = "A",
        threshold = 0.25
      ),
      by_scope = "subject"
    )
  )
  testthat::expect_error(
    annotate_ig_clones(
      idata,
      heavy_threshold = tibble::tibble(
        subject = c("A", "A", "B"),
        threshold = c(0.25, 0.25, 0.25)
      ),
      by_scope = "subject"
    )
  )
  testthat::expect_error(
    annotate_ig_clones(
      idata,
      heavy_threshold = tibble::tibble(
        subject = c("A", "B"),
        threshold = c(0, 0.25)
      ),
      by_scope = "subject"
    )
  )
})


testthat::test_that("heavy mode ignores incompatible light chains", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGL", "TTTT", "IGLV3*01", "IGLJ2*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    mode = "heavy"
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(family_partition(actual, "cell"), "c1,c2")
  testthat::expect_equal(actual$family_id, actual$heavy_family_id)
  testthat::expect_true(all(is.na(actual$light_family_id)))
})


testthat::test_that("paired mode explicitly refines heavy families by light families", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTT", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    mode = "paired"
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell", "heavy_family_id"),
    "c1,c2"
  )
  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1", "c2"))
  )
  testthat::expect_false(anyNA(actual$light_family_id))
})


testthat::test_that("auto mode supports heavy-only receptor schemas", {
  rows <- tibble::tribble(
    ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
    "S1", "c1", "A", "IGH", "AAAA", "IGHV1*01", "IGHJ1*01", 10L,
    "S1", "c2", "A", "IGH", "AAAT", "IGHV1*01", "IGHJ1*01", 10L
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(rows, chains = "IGH"),
    heavy_threshold = 0.25,
    mode = "auto"
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(family_partition(actual, "cell"), "c1,c2")
  testthat::expect_equal(actual$family_id, actual$heavy_family_id)
  testthat::expect_true(all(is.na(actual$light_family_id)))
})


testthat::test_that("heavy mode supports heavy-only receptor schemas", {
  rows <- tibble::tribble(
    ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
    "S1", "c1", "A", "IGH", "AAAA", "IGHV1*01", "IGHJ1*01", 10L,
    "S1", "c2", "A", "IGH", "AAAT", "IGHV1*01", "IGHJ1*01", 10L,
    "S1", "c3", "A", "IGH", "TTTT", "IGHV1*01", "IGHJ1*01", 10L
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(rows, chains = "IGH"),
    heavy_threshold = 0.25,
    mode = "heavy"
  )

  actual <- collect_cell_families(out)

  testthat::expect_equal(
    family_partition(actual, "cell"),
    sort(c("c1,c2", "c3"))
  )
  testthat::expect_equal(actual$family_id, actual$heavy_family_id)
  testthat::expect_true(all(is.na(actual$light_family_id)))
})


testthat::test_that("paired mode rejects heavy-only receptor schemas", {
  rows <- tibble::tribble(
    ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
    "S1", "c1", "A", "IGH", "AAAA", "IGHV1*01", "IGHJ1*01", 10L
  )
  idata <- make_ig_clone_idata(rows, chains = "IGH")

  testthat::expect_error(
    annotate_ig_clones(
      idata,
      heavy_threshold = 0.25,
      mode = "paired"
    ),
    "[Pp]aired mode requires a paired IGH.*IGK.*IGL receptor schema"
  )
})


testthat::test_that("column arguments select nonstandard AIRR columns", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01"
  )
  rows <- paired_ig_rows(cells) |>
    dplyr::rename(
      junction_nt = "junction",
      v_gene = "v_call",
      j_gene = "j_call"
    )

  out <- annotate_ig_clones(
    make_ig_clone_idata(
      rows,
      features = c("junction_nt", "v_gene", "j_gene")
    ),
    heavy_threshold = 0.25,
    junction_col = "junction_nt",
    v_col = "v_gene",
    j_col = "j_gene"
  )

  testthat::expect_equal(
    family_partition(out$annotations, "cell"),
    "c1,c2"
  )
})


testthat::test_that("heavy locus is configurable", {
  rows <- tibble::tribble(
    ~sample, ~cell, ~subject, ~locus, ~junction, ~v_call, ~j_call, ~umi,
    "S1", "c1", "A", "H", "AAAA", "HV1", "HJ1", 10L,
    "S1", "c1", "A", "K", "CCCC", "KV1", "KJ1", 8L,
    "S1", "c2", "A", "H", "AAAT", "HV1", "HJ1", 10L,
    "S1", "c2", "A", "K", "CCCT", "KV1", "KJ1", 8L
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(
      rows,
      chains = c("H", "K|L")
    ),
    heavy_threshold = 0.25,
    heavy_locus = "H"
  )

  testthat::expect_equal(
    family_partition(out$annotations, "cell"),
    "c1,c2"
  )
})


testthat::test_that("family_col renames only the paired family output", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25,
    family_col = "antibody_family"
  )

  testthat::expect_true("antibody_family" %in% colnames(out$annotations))
  testthat::expect_false("family_id" %in% colnames(out$annotations))
  testthat::expect_true(
    all(c("heavy_family_id", "light_family_id") %in%
      colnames(out$annotations))
  )
})


testthat::test_that("rerunning annotation replaces family columns cleanly", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01"
  )
  family_cols <- c(
    "family_id",
    "heavy_family_id",
    "light_family_id"
  )
  idata <- make_ig_clone_idata(paired_ig_rows(cells))

  restrictive <- annotate_ig_clones(
    idata,
    heavy_threshold = 0.20,
    light_threshold = 0.20
  )
  permissive <- annotate_ig_clones(
    restrictive,
    heavy_threshold = 0.25,
    light_threshold = 0.25
  )

  testthat::expect_equal(
    family_partition(restrictive$annotations, "cell"),
    sort(c("c1", "c2"))
  )
  testthat::expect_equal(
    family_partition(permissive$annotations, "cell"),
    "c1,c2"
  )
  testthat::expect_false(any(family_cols %in% colnames(idata$annotations)))
  testthat::expect_true(all(
    vapply(
      family_cols,
      function(column) sum(colnames(permissive$annotations) == column) == 1L,
      logical(1)
    )
  ))
  testthat::expect_false(any(grepl(
    "\\.[xy]$",
    colnames(permissive$annotations)
  )))
})


testthat::test_that("exact receptor identity and observations are preserved", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )
  receptor_col <- immundata::imd_schema("receptor")

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )
  actual <- out$annotations |>
    dplyr::select(
      dplyr::all_of(c(
        "cell",
        "locus",
        receptor_col,
        "family_id",
        "heavy_family_id",
        "light_family_id"
      ))
    ) |>
    dplyr::collect()

  testthat::expect_equal(dplyr::n_distinct(actual[[receptor_col]]), 1L)
  testthat::expect_equal(nrow(actual), 4L)
  testthat::expect_equal(family_partition(actual, "cell"), "c1,c2")
})


testthat::test_that("both chains receive one consistent family annotation", {
  cells <- tibble::tribble(
    ~sample, ~cell, ~subject, ~h_junction, ~h_v_call, ~h_j_call,
    ~l_locus, ~l_junction, ~l_v_call, ~l_j_call,
    "S1", "c1", "A", "AAAA", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01",
    "S1", "c2", "A", "AAAT", "IGHV1*01", "IGHJ1*01",
    "IGK", "CCCT", "IGKV1*01", "IGKJ1*01",
    "S1", "c3", "A", "AATT", "IGHV1*01", "IGHJ1*01",
    "IGK", "TTTT", "IGKV1*01", "IGKJ1*01",
    "S1", "c4", "A", "GGGG", "IGHV3*01", "IGHJ3*01",
    "IGK", "CCCC", "IGKV1*01", "IGKJ1*01"
  )

  out <- annotate_ig_clones(
    make_ig_clone_idata(paired_ig_rows(cells)),
    heavy_threshold = 0.25
  )
  actual <- collect_cell_families(out)
  per_cell <- actual |>
    dplyr::summarise(
      n_loci = dplyr::n_distinct(.data$locus),
      n_family = dplyr::n_distinct(.data$family_id),
      n_heavy_family = dplyr::n_distinct(.data$heavy_family_id),
      n_light_family = dplyr::n_distinct(.data$light_family_id),
      .by = "cell"
    )
  containment <- actual |>
    dplyr::distinct(
      .data$family_id,
      .data$heavy_family_id,
      .data$light_family_id
    ) |>
    dplyr::summarise(
      n_heavy = dplyr::n_distinct(.data$heavy_family_id),
      n_light = dplyr::n_distinct(.data$light_family_id),
      .by = "family_id"
    )
  pair_mapping <- actual |>
    dplyr::distinct(
      .data$family_id,
      .data$heavy_family_id,
      .data$light_family_id
    ) |>
    dplyr::summarise(
      n_family = dplyr::n_distinct(.data$family_id),
      .by = c("heavy_family_id", "light_family_id")
    )

  testthat::expect_true(all(per_cell$n_loci == 2L))
  testthat::expect_true(all(per_cell$n_family == 1L))
  testthat::expect_true(all(per_cell$n_heavy_family == 1L))
  testthat::expect_true(all(per_cell$n_light_family == 1L))
  testthat::expect_true(all(containment$n_heavy == 1L))
  testthat::expect_true(all(containment$n_light == 1L))
  testthat::expect_true(all(pair_mapping$n_family == 1L))
  testthat::expect_s3_class(out$annotations, "duckplyr_df")
})


testthat::test_that("V/J unions precede junction distance and length blocking", {
  for (bridge in c("CCCC", "CCCCC")) {
    rows <- tibble::tibble(
      sample = "S1", cell = c("a", "bridge", "c"), subject = "A",
      locus = "IGH", junction = c("AAAA", bridge, "AAAT"),
      v_call = c("IGHV1*01", " IGHV2*02 | IGHV1*01,IGHV1*02 ", "IGHV2*01"),
      j_call = c("IGHJ1*01", "IGHJ2*01, IGHJ1*02", "IGHJ2*02"),
      umi = 10L
    )
    idata <- make_ig_clone_idata(rows, chains = "IGH")
    before <- idata$annotations |> dplyr::collect()
    out <- annotate_ig_clones(idata, heavy_threshold = 0.25, mode = "heavy")
    testthat::expect_equal(
      family_partition(out$annotations, "cell"), c("a,c", "bridge")
    )
    testthat::expect_equal(
      out$annotations |> dplyr::select(dplyr::all_of(names(before))) |>
        dplyr::arrange(cell) |> dplyr::collect(),
      before |> dplyr::arrange(cell)
    )
    testthat::expect_s3_class(out$annotations, "duckplyr_df")
  }
})


testthat::test_that("ambiguity bridges stay within scopes and require both genes", {
  rows <- tibble::tibble(
    sample = "S1", cell = c("a", "c", "foreign", "wrong_j"),
    subject = c("A", "A", "B", "A"), locus = "IGH",
    junction = "AAAA",
    v_call = c("IGHV1", "IGHV2", "IGHV1,IGHV2", "IGHV1,IGHV2"),
    j_call = c("IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2"), umi = 10L
  )
  out <- annotate_ig_clones(
    make_ig_clone_idata(rows, chains = "IGH"), heavy_threshold = 0.25
  )
  testthat::expect_equal(
    family_partition(out$annotations, "cell"), sort(rows$cell)
  )
})


testthat::test_that("equivalent calls normalize and empty calls cannot connect", {
  rows <- tibble::tibble(
    sample = "S1", cell = c("a", "b", "empty1", "empty2"), subject = "A",
    locus = "IGH", junction = "AAAA",
    v_call = c("IGHV1*01", " IGHV1*02 | IGHV1*01 ", ",|", ",|"),
    j_call = "IGHJ1", umi = 10L
  )
  # Different junctions ensure the incomplete calls represent different receptors.
  rows$junction[4] <- "AAAT"
  out <- annotate_ig_clones(
    make_ig_clone_idata(rows, chains = "IGH"), heavy_threshold = 0.25
  )
  testthat::expect_equal(
    family_partition(out$annotations, "cell"), c("a,b", "empty1", "empty2")
  )
})


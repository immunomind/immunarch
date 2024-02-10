if (is.environment(immdata)) {
  test_data = immdata
} else {
  test_data = immdata$data
}

test_that("snapshot_repExplore_count", {
  local_edition(3)
  expect_snapshot({
    repExplore(test_data, .method = "count", .col = "nt", .coding = TRUE) %>%
      head()
  })
  expect_snapshot({
    repExplore(test_data, .method = "count", .col = "nt", .coding = TRUE) %>%
      tail()
  })
})

test_that("snapshot_repExplore_clones", {
  local_edition(3)
  expect_snapshot({
    repExplore(immdata$data, .method = "clones", .col = "nt", .coding = TRUE) %>%
      head()
  })
  expect_snapshot({
    repExplore(immdata$data, .method = "clones", .col = "nt", .coding = TRUE) %>%
      tail()
  })
})

test_that("snapshot_repClonality_homeo", {
  local_edition(3)
  expect_snapshot({
    repClonality(test_data, .method = "homeo")
  })
})

test_that("snapshot_repDiversity_chao1", {
  local_edition(3)
  expect_snapshot({
    repDiversity(test_data, .method = "chao1", .verbose = FALSE)
  })
})

test_that("snapshot_repDiversity_D50", {
  local_edition(3)
  expect_snapshot({
    repDiversity(test_data, .method = "D50", .verbose = FALSE)
  })
})

test_that("snapshot_repDiversity_raref", {
  local_edition(3)
  expect_snapshot({
    repDiversity(test_data, .method = "raref", .verbose = FALSE) %>%
      group_by(Sample) %>%
      slice_tail(n=2) %>%
      ungroup
  })
})

test_that("snapshot_geneUsage_segment", {
  local_edition(3)
  expect_snapshot({
    geneUsage(test_data, .gene = "hs.trbj", .type = "segment")
  })
})

test_that("snapshot_geneUsage_family", {
  local_edition(3)
  expect_snapshot({
    geneUsage(test_data, .gene = "hs.trbj", .type = "family")
  })
})

test_that("snapshot_repOverlap_shared", {
  local_edition(3)
  expect_snapshot({
    repOverlap(test_data, .method = "shared")
  })
})

test_that("snapshot_repOverlap_jaccard", {
  local_edition(3)
  expect_snapshot({
    repOverlap(test_data, .method = "jaccard", .verbose = FALSE)
  })
})

test_that("snapshot_repOverlap_morisita", {
  local_edition(3)
  expect_snapshot({
    repOverlap(test_data, .method = "morisita", .verbose = FALSE)
  })
})

test_that("snapshot_publicRepertoire", {
  local_edition(3)
  expect_snapshot({
    publicRepertoire(test_data, .min.samples = 5, .verbose = FALSE)
  })
})

test_that("snapshot_trackClonotypes", {
  local_edition(3)
  expect_snapshot({
    trackClonotypes(
      test_data,
      c(
        "CASSLEETQYF",
        "CASSFQETQYF",
        "CASSDSSGGANEQFF",
        "CASSLGETQYF"
      ),
      .col = "aa"
    )
  })
})

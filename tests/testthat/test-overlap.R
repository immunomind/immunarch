for (column in c("nt", "aa", "nt+v", "aa+v", "aa+v+j")) {
  for (method in c("public", "overlap", "jaccard", "tversky", "cosine", "morisita", "inc+public", "inc+jaccard", "inc+morisita")) {
    test_name <- paste0("method:", method, ":", column)

    compute_res <- apply_DF_DT(frame_data, table_data,
      repOverlap,
      .method = method, .col = column, .verbose = F, .verbose.inc = F
    )
    test_that(ap(test_name, "compute_"), {
      expect_equal(compute_res[[1]], compute_res[[2]])
    })

    # vis_res = vis_results(compute_res)
    # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
  }
}

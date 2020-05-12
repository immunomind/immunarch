for (column in c("nt", "aa", "nt+v", "aa+v", "nt+j", "aa+j")) {
  for (quant in c("id", "count")) {
    test_name <- paste0("quant:", quant, ":", column)

    compute_res <- apply_DF_DT(frame_data[[1]], table_data[[1]],
      spectratype,
      .quant = quant, .col = column
    )
    test_that(ap(test_name, "compute_"), {
      expect_equal(compute_res[[1]], compute_res[[2]])
    })

    # vis_res = vis_results(compute_res)
    # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
  }
}

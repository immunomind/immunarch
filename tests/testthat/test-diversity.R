for (column in c("nt", "aa", "nt+v", "aa+v", "aa+v+j")) {
  for (method in c("chao1", "hill", "div", "gini.simp", "inv.simp", "gini", "raref", "dxx", "d50")) {
    if (method == "raref") {
      for (norm_val in c(TRUE, FALSE)) {
        test_name <- paste0("method:", method, ":", column, ".norm:", norm_val)

        compute_res <- apply_DF_DT(frame_data, table_data,
          repDiversity,
          .method = method, .col = column, .verbose = F, .norm = norm_val
        )
        test_that(ap(test_name, "compute_"), {
          expect_equal(compute_res[[1]], compute_res[[2]])
        })

        # vis_res = vis_results(compute_res)
        # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
      }
    } else if (method == "dxx") {
      for (perc_val in c(10, 25, 75)) {
        test_name <- paste0("method:", method, ":", column, ".perc:", perc_val)

        compute_res <- apply_DF_DT(frame_data, table_data,
          repDiversity,
          .method = method, .col = column, .verbose = F, .perc = perc_val
        )
        test_that(ap(test_name, "compute_"), {
          expect_equal(compute_res[[1]], compute_res[[2]])
        })

        # vis_res = vis_results(compute_res)
        # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
      }
    } else {
      test_name <- paste0("method:", method, ":", column)

      compute_res <- apply_DF_DT(frame_data, table_data,
        repDiversity,
        .method = method, .col = column, .verbose = F
      )
      test_that(ap(test_name, "compute_"), {
        expect_equal(compute_res[[1]], compute_res[[2]])
      })

      # vis_res = vis_results(compute_res)
      # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
    }
  }
}

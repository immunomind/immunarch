for (column in c("nt", "aa", "nt+v", "aa+v", "aa+v+j")) {
  for (method in c("volume", "count", "len", "clones")) {
    for (coding_val in c(TRUE, FALSE)) {
      if (method == "len") {
        if (column %in% c("nt", "aa")) {
          test_name <- paste0("method:", method, ":", column, ".coding:", coding_val)

          compute_res <- apply_DF_DT(frame_data, table_data,
            repExplore,
            .method = method, .col = column, .coding = coding_val
          )
          test_that(ap(test_name, "compute_"), {
            expect_equal(compute_res[[1]], compute_res[[2]])
          })

          # vis_res = vis_results(compute_res)
          # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
        }
      } else {
        test_name <- paste0("method:", method, ":", column, ".coding:", coding_val)

        compute_res <- apply_DF_DT(frame_data, table_data,
          repExplore,
          .method = method, .col = column, .coding = coding_val
        )
        test_that(ap(test_name, "compute_"), {
          expect_equal(compute_res[[1]], compute_res[[2]])
        })

        # vis_res = vis_results(compute_res)
        # test_that(ap(test_name, "vis_"), { expect_equal(vis_res[[1]], vis_res[[2]]) })
      }
    }
  }
}

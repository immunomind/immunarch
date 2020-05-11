for (quant in c(NA, "count")) {
  for (type in c("segment", "allele", "family")) {
    for (ambig in c("inc", "exc", "maj")) {
      for (norm_val in c(TRUE, FALSE)) {
        test_name <- paste0("quant:", quant, "_type:", type, "_ambig:", ambig, "_norm:", norm_val)

        compute_res <- apply_DF_DT(frame_data, table_data,
          geneUsage,
          .quant = quant, .ambig = ambig, .type = type, .norm = norm_val
        )
        test_that(ap(test_name, "compute_"), {
          expect_equal(compute_res[[1]], compute_res[[2]])
        })
      }
    }
  }
}

skip_on_ci()
skip_on_bioc()
skip_if_offline("rest.ensembl.org")

test_that("single and multiple getGene() calls are identical", {
  expect_warning(
    multi_call <- getGene(
      id = c("ENSG00000157764", "ENSG00000004939"),
      type = "ensembl_gene_id"
    ),
    "'band' column information is not available"
  )

  expect_warning(
    single_call_1 <- getGene(
      id = "ENSG00000157764",
      type = "ensembl_gene_id"
    ),
    "'band' column information is not available"
  )

  expect_warning(
    single_call_2 <- getGene(
      id = "ENSG00000004939",
      type = "ensembl_gene_id"
    ),
    "'band' column information is not available"
  )

  expect_identical(
    multi_call,
    rbind(single_call_1, single_call_2)
  )

  expect_named(
    multi_call,
    c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "description",
      "chromosome_name",
      "band",
      "strand",
      "start_position",
      "end_position"
    )
  )
})

test_that("getGene() fails with unsupported ID types", {
  expect_error(
    getGene(
      id = "ENST00000357654",
      type = "unknown_id"
    ),
    "Only Ensembl Gene IDs"
  )
})

skip_on_ci()
skip_on_bioc()
skip_if_offline("rest.ensembl.org")

test_that("getLDS() works with homologies", {
  res <- getLDS(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = c("ENSG00000157764", "ENSG00000004939"),
    attributesL = c("ensembl_gene_id"),
    speciesL = "mouse"
  ) |>
    expect_no_error() |>
    expect_no_warning()

  expect_named(
    res,
    c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "external_gene_name",
      "ensembl_gene_id.1"
    )
  )
})

test_that("getLDS() works with gene symbols", {
  res <- getLDS(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
    filters = "external_gene_name",
    values = c("APOE", "MAPT"),
    species = "human",
    attributesL = c("ensembl_gene_id"),
    speciesL = "mouse"
  ) |>
    expect_no_error() |>
    expect_no_warning()

  expect_named(
    res,
    c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "external_gene_name",
      "ensembl_gene_id.1"
    )
  )
})

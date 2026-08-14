#' List supported attributes
#'
#' @param ... Ignored. Used to catch no longer necessary parameters from the
#'   \pkg{biomaRt} functions.
#'
#' @returns A character vector of supported attributes to be used in [getBM()]
#'   and [getLDS()].
#'
#' @export
#'
#' @examples
#' listAttributes()
#'
listAttributes <- function(...) {
  c(
    .listGeneLevelAttributes(),
    .listTranscriptLevelAttributes()
  ) |>
    unique() |>
    sort(method = "radix")
}

.listGeneLevelAttributes <- function() {
  c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "transcript_biotype",
    "version",
    "hgnc_symbol"
  )
}

.listTranscriptLevelAttributes <- function() {
  c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_peptide_id",
    "transcript_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "version"
  )
}

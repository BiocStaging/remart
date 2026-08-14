#' List supported attributes
#'
#' @param ... Ignored. Used to catch no longer necessary parameters from the
#'   \pkg{biomaRt} functions.
#'
#' @export
#'
#' @examples
#' listAttributes()
#'
listAttributes <- function(...) {
  transcript_level_attributes <- c(
    "ensembl_transcript_id",
    "ensembl_peptide_id",
    "transcript_biotype"
  )
  supported_attributes <- c(
    transcript_level_attributes,
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
  return(supported_attributes)
}

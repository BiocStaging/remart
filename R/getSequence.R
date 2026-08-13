#' Get Sequences from Ensembl
#' 
#' @inheritParams biomaRt::getSequence
#' @inheritParams getGene
#' 
#' @export
#' 
#' @examples
#' remart::getSequence(
#'   seqType = "gene_exon_intron",
#'   type = "ensembl_gene_id",
#'   id = "ENSG00000001497"
#' )
#'
#' ids <- c(
#'  "ENSG00000003987",
#'  "ENSG00000004939"
#' )
#' remart::getSequence(
#'   seqType = "gene_exon_intron",
#'   type = "ensembl_gene_id",
#'   id = ids
#' )
getSequence <- function(
  chromosome,
  start,
  end,
  id,
  type,
  seqType,
  upstream,
  downstream,
  ...
) {
  if (!missing(chromosome)) {
    stop("Not implemented yet")
  }

  if (!missing(id)) {
    if (type %notin% c("ensembl_gene_id", "ensembl_transcript_id")) {
      stop("Only Ensembl IDs (ENS...) are supported at the moment")
    }
    if (!is.list(id) && length(id) == 1) {
      id <- list(id)
    }

    req <- list(
      ids = id
    )
    if (!missing(upstream)) {
      req$expand_5prime <- upstream
    }
    if (!missing(downstream)) {
      req$expand_3prime <- downstream
    }

    req$type <- switch(
      seqType,
      "gene_exon_intron" = "genomic",
      "cdna" = "cdna",
      "peptide" = "protein", # see example use in DominoEffect package
      stop("Invalid seqType. Must be one of 'gene_exon_intron', 'cdna', or 'peptide'.")
    )

    httr2::request("https://rest.ensembl.org/") |> 
      httr2::req_url_path("sequence/id") |> 
      httr2::req_method("POST") |> 
      httr2::req_user_agent(REMART_USER_AGENT) |> 
      httr2::req_body_json(req) |> 
      httr2::req_perform() |>
      httr2::resp_body_json(simplifyVector = TRUE)
  }
}
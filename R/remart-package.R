#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Backport from R 4.6.0
`%notin%` <- function(x, table) {
  match(x, table, nomatch = 0L) == 0
}

# Backport of the base R null-coalescing operator (R >= 4.4.0)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Batch-fetch Ensembl feature metadata via the REST `lookup/id` endpoint
#'
#' @param ids character vector of Ensembl stable IDs (genes, transcripts or
#'   translations are all supported by this endpoint).
#' @param expand if `TRUE`, also return the child features (e.g. transcripts
#'   and translations for a gene).
#'
#' @return A named list of parsed JSON objects, keyed by `ids`. Missing IDs
#'   are returned as `NULL` entries.
#' @noRd
.remart_lookup_id <- function(ids, expand = FALSE) {
  ids <- unique(as.character(ids))

  req <- httr2::request("https://rest.ensembl.org") |>
    httr2::req_url_path("/lookup/id") |>
    httr2::req_method("POST") |>
    httr2::req_user_agent("remart R package") |>
    httr2::req_body_json(list(ids = as.list(ids)))

  if (expand) {
    req <- httr2::req_url_query(req, expand = 1)
  }

  res <- req |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = FALSE)

  # Ensure every requested id has an (possibly NULL) entry, and that the
  # result is returned in the same order as `ids`.
  res[ids]
}
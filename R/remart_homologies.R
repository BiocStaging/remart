#' Fetch orthologues of `id` in `target_species` from the REST `homology/id`
#' endpoint
#' @noRd
.remart_homologies <- function(id, species, target_species) {
  resp <- httr2::request("https://rest.ensembl.org") |>
    httr2::req_url_path(sprintf("/homology/id/%s/%s", species, id)) |>
    httr2::req_url_query(
      type = "orthologues",
      target_species = target_species,
      `content-type` = "application/json"
    ) |>
    httr2::req_user_agent(REMART_USER_AGENT) |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = FALSE)

  resp$data[[1]]$homologies
}

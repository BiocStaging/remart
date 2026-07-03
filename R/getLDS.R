#' Retrieve orthologous gene annotations across two species
#'
#' @inheritParams biomaRt::getLDS
#' @param speciesL Ensembl name (e.g. `"mouse"` or `"mus_musculus"`) of the
#'   species to look up orthologues in. This replaces the `martL` argument
#'   used in \pkg{biomaRt}, as connections to Ensembl datasets are no longer
#'   needed, but the species of the linked dataset can not be inferred from
#'   `attributesL`/`filtersL`/`valuesL` alone.
#' @param ... Ignored. Used to catch no longer necessary parameters such as
#'   `mart`, `martL`, `verbose`, `uniqueRows` and `bmHeader` from
#'   \pkg{biomaRt} functions.
#'
#' @details
#' This function relies on the Ensembl REST `homology/id` endpoint, and
#' therefore only supports retrieving orthologues, and not the more general
#' cross-database attribute linking `biomaRt::getLDS()` could perform.
#'
#' Supported filters/attributes (for both `filters`/`attributes` and
#' `filtersL`/`attributesL`): `ensembl_gene_id`, `ensembl_peptide_id`,
#' `external_gene_name`, `description`, `chromosome_name`,
#' `start_position`, `end_position`, `strand`, `gene_biotype`.
#'
#' Only `filters = "ensembl_gene_id"` is supported.
#'
#' @export
#'
#' @examples
#' remart::getLDS(
#'   attributes = c("ensembl_gene_id", "external_gene_name"),
#'   filters = "ensembl_gene_id",
#'   values = "ENSG00000157764",
#'   attributesL = c("ensembl_gene_id", "external_gene_name"),
#'   speciesL = "mouse"
#' )
getLDS <- function(
  attributes,
  filters = "",
  values = "",
  attributesL,
  filtersL = "",
  valuesL = "",
  speciesL,
  ...
) {
  supported_attributes <- c(
    "ensembl_gene_id", "ensembl_peptide_id", "external_gene_name",
    "description", "chromosome_name", "start_position", "end_position",
    "strand", "gene_biotype"
  )

  if (length(filters) != 1 || filters != "ensembl_gene_id") {
    stop("Only filters = \"ensembl_gene_id\" is supported at the moment.")
  }

  if (length(filtersL) != 1 || filtersL %notin% c("", supported_attributes)) {
    stop(
      "filtersL must be \"\" or one of: ",
      paste(supported_attributes, collapse = ", ")
    )
  }

  unsupported_attributes <- setdiff(c(attributes, attributesL), supported_attributes)
  if (length(unsupported_attributes) > 0) {
    stop(
      "Unsupported attribute(s): ",
      paste(unsupported_attributes, collapse = ", "),
      ". Supported attributes are: ",
      paste(supported_attributes, collapse = ", ")
    )
  }

  if (missing(speciesL) || length(speciesL) != 1 || speciesL == "") {
    stop(
      "`speciesL` (e.g. \"mouse\" or \"mus_musculus\"), identifying the ",
      "target species, must be provided. It replaces the `martL` argument ",
      "used in biomaRt."
    )
  }

  values <- as.character(values)
  values <- values[values != ""]
  if (length(values) == 0) {
    stop("`values` must contain at least one identifier.")
  }

  source_genes <- .remart_lookup_id(values, expand = FALSE)

  missing_ids <- values[vapply(source_genes, is.null, logical(1))]
  if (length(missing_ids) > 0) {
    warning(
      "The following identifiers were not found and will be ignored: ",
      paste(missing_ids, collapse = ", ")
    )
  }

  rows <- lapply(values, function(id) {
    gene <- source_genes[[id]]
    if (is.null(gene)) {
      return(NULL)
    }

    homologies <- .remart_homologies(
      id,
      species = gene$species,
      target_species = speciesL
    )
    if (length(homologies) == 0) {
      return(NULL)
    }

    target_ids <- vapply(homologies, function(h) h$target$id, character(1))
    target_genes <- .remart_lookup_id(target_ids, expand = FALSE)

    lapply(homologies, function(h) {
      target_gene <- target_genes[[h$target$id]]

      if (filtersL != "") {
        filter_value <- .remart_bm_attr(filtersL, gene = target_gene)
        if (filter_value %notin% valuesL) {
          return(NULL)
        }
      }

      row_source <- .remart_bm_row(attributes, gene = gene)
      if ("ensembl_peptide_id" %in% attributes) {
        row_source$ensembl_peptide_id <- h$source$protein_id %||% NA_character_
      }

      row_target <- .remart_bm_row(attributesL, gene = target_gene)
      if ("ensembl_peptide_id" %in% attributesL) {
        row_target$ensembl_peptide_id <- h$target$protein_id %||% NA_character_
      }
      names(row_target) <- paste0(names(row_target), ".1")

      cbind(row_source, row_target)
    })
  })

  rows <- Filter(Negate(is.null), unlist(rows, recursive = FALSE))

  cols <- c(attributes, paste0(attributesL, ".1"))
  if (length(rows) == 0) {
    empty_cols <- stats::setNames(
      replicate(length(cols), character(0), simplify = FALSE),
      cols
    )
    return(as.data.frame(empty_cols))
  }

  df <- unique(do.call(rbind, rows))
  rownames(df) <- NULL
  df
}

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
    httr2::req_user_agent("remart R package") |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = FALSE)

  resp$data[[1]]$homologies
}
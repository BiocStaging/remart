#' Retrieve gene/transcript annotations from Ensembl
#'
#' @inheritParams biomaRt::getBM
#' @param ... Ignored. Used to catch no longer necessary parameters such as
#'   `mart`, `checkFilters`, `verbose`, `uniqueRows`, `bmHeader`, `quote` and
#'   `useCache` from \pkg{biomaRt} functions.
#'
#' @details
#' Only a subset of the attributes and filters supported by the `biomaRt`
#' package are currently implemented, as this data has to be retrieved
#' through the Ensembl REST `lookup/id` endpoint rather than a generic
#' BioMart query engine.
#'
#' Supported filters (only one can be used at a time): `ensembl_gene_id`,
#' `ensembl_transcript_id`.
#'
#' Supported attributes: `ensembl_gene_id`, `ensembl_transcript_id`,
#' `ensembl_peptide_id`, `external_gene_name`, `description`,
#' `chromosome_name`, `start_position`, `end_position`, `strand`,
#' `gene_biotype`, `transcript_biotype`, `version`.
#'
#' @export
#'
#' @examples
#' remart::getBM(
#'   attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
#'   filters = "ensembl_gene_id",
#'   values = c("ENSG00000157764", "ENSG00000004939")
#' )
getBM <- function(
  attributes,
  filters = "",
  values = "",
  ...
) {
  stopifnot(
    is.character(attributes),
    is.character(filters),
    is.character(values)
  )

  supported_filters <- c("ensembl_gene_id", "ensembl_transcript_id")
  transcript_level_attributes <- c(
    "ensembl_transcript_id", "ensembl_peptide_id", "transcript_biotype"
  )
  supported_attributes <- c(
    transcript_level_attributes,
    "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id",
    "external_gene_name", "description", "chromosome_name",
    "start_position", "end_position", "strand",
    "gene_biotype", "transcript_biotype", "version"
  )

  if (length(filters) != 1L || filters %notin% supported_filters) {
    stop(
      "Only a single filter is supported at the moment, and must be one of: ",
      toString(supported_filters)
    )
  }

  unsupported_attributes <- setdiff(attributes, supported_attributes)
  if (length(unsupported_attributes) > 0L) {
    stop(
      "Unsupported attribute(s): ",
      toString(unsupported_attributes),
      ".\nSupported attributes are: ",
      toString(supported_attributes),
      ".\nPlease file an issue at ",
      "https://github.com/Huber-group-EMBL/remart/issues",
      " if you would like us to add support for new attributes."
    )
  }

  values <- values[nzchar(values)]
  if (length(values) == 0L) {
    stop("`values` must contain at least one identifier.")
  }

  needs_transcripts <- any(attributes %in% transcript_level_attributes)

  if (filters == "ensembl_gene_id") {
    genes <- .remart_lookup_id(values, expand = needs_transcripts)

    missing_ids <- values[lengths(genes) == 0L]
    if (length(missing_ids) > 0L) {
      warning(
        "The following identifiers were not found and will be ignored: ",
        toString(missing_ids)
      )
    }

    rows <- lapply(genes, function(gene) {
      if (is.null(gene)) {
        return(NULL)
      }
      transcripts <- gene$Transcript
      if (!needs_transcripts || is.null(transcripts)) {
        list(.remart_bm_row(attributes, gene = gene))
      } else {
        lapply(
          transcripts,
          function(x) .remart_bm_row(attributes, gene = gene, transcript = x)
        )
      }
    })
  } else {
    transcripts <- .remart_lookup_id(values, expand = TRUE)

    missing_ids <- values[vapply(transcripts, is.null, logical(1))]
    if (length(missing_ids) > 0) {
      warning(
        "The following identifiers were not found and will be ignored: ",
        toString(missing_ids)
      )
    }

    needs_genes <- any(attributes %notin% transcript_level_attributes)
    gene_ids <- unique(unlist(lapply(transcripts, `[[`, "Parent")))
    genes <- if (needs_genes) .remart_lookup_id(gene_ids, expand = FALSE) else list()

    rows <- lapply(transcripts, function(transcript) {
      if (is.null(transcript)) {
        return(NULL)
      }
      gene <- if (needs_genes) genes[[transcript$Parent]] else NULL
      list(.remart_bm_row(attributes, gene = gene, transcript = transcript))
    })
  }

  rows <- unlist(rows, recursive = FALSE)

  if (length(rows) == 0) {
    empty_cols <- stats::setNames(
      replicate(length(attributes), character(0), simplify = FALSE),
      attributes
    )
    return(as.data.frame(empty_cols))
  }

  df <- unique(do.call(rbind, rows))
  rownames(df) <- NULL
  df
}

#' Build a single-row data.frame of `attributes` from parsed lookup/id objects
#' @noRd
.remart_bm_row <- function(attributes, gene = NULL, transcript = NULL) {
  values <- lapply(attributes, .remart_bm_attr, gene = gene, transcript = transcript)
  names(values) <- attributes
  as.data.frame(values, stringsAsFactors = FALSE)
}

#' Extract the value of a single `biomaRt`-style attribute from parsed
#' `lookup/id` objects
#' @noRd
.remart_bm_attr <- function(attr, gene = NULL, transcript = NULL) {
  translation <- transcript$Translation
  main <- gene %||% transcript

  switch(
    attr,
    ensembl_gene_id = gene$id %||% transcript$Parent %||% NA_character_,
    ensembl_transcript_id = transcript$id %||% NA_character_,
    ensembl_peptide_id = translation$id %||% NA_character_,
    external_gene_name = gene$display_name %||% NA_character_,
    description = gene$description %||% NA_character_,
    chromosome_name = main$seq_region_name %||% NA_character_,
    start_position = main$start %||% NA_integer_,
    end_position = main$end %||% NA_integer_,
    strand = main$strand %||% NA_integer_,
    gene_biotype = gene$biotype %||% NA_character_,
    transcript_biotype = transcript$biotype %||% NA_character_,
    version = (transcript %||% gene)$version %||% NA_integer_,
    stop("Unsupported attribute: ", attr)
  )
}
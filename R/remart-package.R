#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats setNames
## usethis namespace: end
NULL

# Backport from R 4.4.0
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Backport from R 4.6.0
`%notin%` <- function(x, table) {
  match(x, table, nomatch = 0L) == 0
}

#' Check if input is a numeric scalar, vector, or matrix, with no missing values
#'
#' @description
#' This function checks to see if `x` is a numeric object with non-missing or
#' non-infinite entries. `x` may be a scalar, a vector, or a matrix. If `x` is
#' one of these or has any missing or infinite (if specified) or non-numeric
#' entries, the function stops the current R process and returns an error.
#'
#' @param x Any R object.
#' @param object_type A character object, one of `scalar`, `vector`, or `matrix` (default `scalar`), indicating which data value type is expected. Note, for this function, a `1`\eqn{\times}`1` matrix is not a scalar, and an `n`\eqn{\times}`1` matrix is not a vector.
#' @param finite A logical `TRUE` or `FALSE` (default `TRUE`), indicating whether the function should check for `-Inf` or `Inf`.
#'
#' @return No returned value
#'
check_numeric <- function(x, object_type = c("scalar", "vector", "matrix")[1], finite = TRUE) {
  # Get variable name:
  name_x <- deparse(substitute(x))

  # Determine whether "stop" message included finite entries comment:
  finite_label <- if (finite) " with finite entries" else ""

  # Determine "stop" message:
  stop_message <- sprintf("%s must be a numeric %s%s.", name_x, object_type, finite_label)


  # Check for numeric data type:
  bool_check <- tryCatch(mode(x) != "numeric", error = function(e) TRUE, warning = function(w) TRUE)
  if (bool_check) stop(stop_message)

  # Check for vector or matrix structure:
  if (object_type != "matrix") {
    bool_check <- tryCatch(!is.vector(x), error = function(e) TRUE, warning = function(w) TRUE)
  } else {
    bool_check <- tryCatch(!is.matrix(x), error = function(e) TRUE, warning = function(w) TRUE)
  }
  if (bool_check) stop(stop_message)

  # If scalar, check for length = 1:
  if (object_type == "scalar") {
    bool_check <- tryCatch(length(x) != 1, error = function(e) TRUE, warning = function(w) TRUE)
    if (bool_check) stop(stop_message)
  }

  # If finite, check for finite (else just check for non-NA/NaN):
  if (finite) {
    bool_check <- tryCatch(!any(is.finite(x)), error = function(e) TRUE, warning = function(w) TRUE)
    if (bool_check) stop(stop_message)
  } else {
    bool_check <- tryCatch(any(is.na(x)), error = function(e) TRUE, warning = function(w) TRUE)
    if (bool_check) stop(stop_message)
  }
}

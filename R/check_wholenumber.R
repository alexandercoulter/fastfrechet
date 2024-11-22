#' Check if input is a whole number (positive integer)
#'
#' @param x 
#'
#' @return No return output.
#' 
check_wholenumber = function(x){
  
  # Get variable name:
  name_x = deparse(substitute(x))
  
  # Determine "stop" message:
  stop_message = sprintf("%s must be a positive integer.", name_x)
  
  # Initial bool_check FALSE:
  bool_check = FALSE
  
  
  # Check for vector:
  bool_check = tryCatch(!is.vector(x), error = function(e) TRUE, warning = function(w) TRUE)
  if(bool_check) stop(stop_message)
  
  # Check for length = 1:
  bool_check = tryCatch(length(x) != 1, error = function(e) TRUE, warning = function(w) TRUE)
  if(bool_check) stop(stop_message)
  
  # Check for finiteness:
  bool_check = tryCatch(!is.finite(x), error = function(e) TRUE, warning = function(w) TRUE)
  if(bool_check) stop(stop_message)
  
  # Check for integer:
  bool_check = tryCatch(x != as.integer(x), error = function(e) TRUE, warning = function(w) TRUE)
  if(bool_check) stop(stop_message)
  
  # Check for > 0:
  bool_check = tryCatch(x < 1, error = function(e) TRUE, warning = function(w) TRUE)
  if(bool_check) stop(stop_message)
  
}

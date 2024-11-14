#' Title
#'
#' @param x 
#' @param finite 
#' @param object_type 
#'
#' @return
#'
#' @examples
check_numeric = function(x, object_type = c("scalar", "vector", "matrix"), finite = TRUE){
  
  # Get variable name:
  name_x = deparse(substitute(x))
  
  # Determine whether "stop" message included finite entries comment:
  finite_label = if(finite) " with finite entries" else ""
  
  # Determine "stop" message:
  stop_message = sprintf("%s must be a numeric %s%s.", name_x, finite_label, object_type)
  
  
  # Check for numeric data type:
  bool_check = tryCatch(mode(x) != "numeric", error = function(e) TRUE)
  if(bool_check) stop(stop_message)
  
  # Check for vector or matrix structure:
  if(object_type != "matrix"){
    
    bool_check = tryCatch(!is.vector(x), error = function(e) TRUE)
    
  } else {
    
    bool_check = tryCatch(!is.matrix(x), error = function(e) TRUE)
    
  }
  if(bool_check) stop(stop_message)
  
  # If scalar, check for length = 1:
  if(scalar){
    
    bool_check = tryCatch(length(x) != 1, error = function(e) TRUE)
    if(bool_check) stop(stop_message)
    
  }
  
  # If finite, check for finite (else just check for non-NA/NaN):
  if(finite){
    
    bool_check = tryCatch(!any(is.finite(x)), error = function(e) TRUE)
    if(bool_check) stop(stop_message)
    
  } else {
    
    bool_check = tryCatch(is.na(x), error = function(e) TRUE)
    if(bool_check) stop(stop_message)
    
  }
  
}

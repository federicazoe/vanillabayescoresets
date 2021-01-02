#' Helper function for printing status messages
#'
#' @param msg a character containing the message to print.
#' @param verbose a boolean for whether or not the msg should be printed.
#'
#' @return NA.
#' 
#' @noRd
#'
#' @examples
#' print_verbose("hello world", TRUE)
print_verbose <- function(msg, verbose) {
  if(verbose) {
    print(msg)
  }
}
#' Combine a list of vectors
#'
#' Intended to combine vectors where, for each position, only one of the vectors contains data (i.e. the remaining are NA's).
#' 
#' Copied from actel
#' 
#' @param input a list of vectors with non-overlapping data.
#' 
#' @return A single vector where all data has been combined.
#' 
#' @keywords internal
#' 
combine <- function(input) {
  if (!inherits(input, "list")) 
    stop("'combine' is only intended to combine a list of vectors to a single vector.")
  if (length(input) == 1) {
    output <- input[[1]]
  } else {
    if (var(unlist(lapply(input, length))) != 0) 
      stop("All vectors to combine should have the same length.")
    output <- input[[1]]
    for (i in 2:length(input)) {
      to.replace <- !is.na(input[[i]])
      if (any(!is.na(output)[to.replace])) 
        stop("Trying to combine value to an already used position.")
      output[to.replace] <- input[[i]][to.replace]
    }
  }
  return(output)
}

#' Remove Code Spaces from transmitter names
#' 
#' Copied from actel
#' 
#' @param input A vector of transmitter names
#' 
#' @keywords internal
#' 
stripCodeSpaces <- function(input) {
  unlist(lapply(input, function(x) tail(unlist(strsplit(x, "-")), 1)))
}

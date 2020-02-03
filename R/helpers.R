#' Convert numeric time to HH:MM
#' 
#' Copied from \link[actel]{actel}
#'
#' @param x Single string or a vector of strings containing hours:minutes or hours:minutes:seconds.
#' @param format the format of x, one of "h" (hours), "m", (minutes) or "s" (seconds).
#' @param seconds Logical; If TRUE, output is returned in HH:MM:SS format.
#' 
#' @return Decimal hour equivalent (single value or vector)
#' 
#' @keywords internal
#' 
minuteTime <- function(x, format = c("h", "m", "s"), seconds = TRUE) {
  format <- match.arg(format)
  .converter <- function(x) {
    if(!is.na(x)){
      if(x < 0){
        x <- abs(x)
        neg = TRUE
      } else neg = FALSE
      if(format == "h") 
        x = x
      if(format == "m") 
        x = x/60
      if(format == "s") 
        x = x/3600
      m = x %% 1
      h = x - m
      m = 60 * m
      s = m %% 1
      m = m - s
      s = round(60 * s, 0)
      if (h < 10) h <- paste0(0, h)
      if (!seconds & s>30) m = m + 1
      if (m < 10) m <- paste0(0, m)
      if (s < 10) s <- paste0(0, s)
      if (seconds) 
        x <- paste(h, m, s, sep = ":")
      else 
        x <- paste(h, m, sep = ":")
      if (neg) x <- paste("-", x)
    }
    return(x)
  }
  if (length(x) < 1) stop("Input appears to be empty.")
  if (!is.numeric(x)) stop("Input is not numeric.")
  if (length(x) == 1) output <- .converter(x)
  if (length(x) > 1) output <- unlist(lapply(x, .converter))
  return(output)
}


#' Forcefully round a number up
#'
#' Forces the rounding of the input to the next higher rounded value.
#' 
#' Copied from \link[actel]{actel}
#' 
#' @param input The value to be rounded.
#' @param to The level of rounding to be applied (i.e. to=10 will round 14.2 to 20; to=1 will round i to 15).
#' 
#' @return The rounded value
#' 
#' @keywords internal
#' 
roundUp <- function(input, to = 10) {
  if (inherits(input, "list"))
    lapply(input, function(input) to * (input %/% to + as.logical(input %% to)))
  else
    to * (input %/% to + as.logical(input %% to))
}

#' Forcefully round a number down
#'
#' Forces the rounding of the input to the next lower rounded value.
#' 
#' @param input The value to be rounded.
#' @param to The level of rounding to be applied (i.e. to=10 will round 14.8 to 10; to=1 will round i to 14).
#' 
#' @return The rounded value
#' 
#' @keywords internal
#' 
roundDown <- function(input, to = 10) {
  to * (input%/%to)
}


#' Combine a list of vectors
#'
#' Intended to combine vectors where, for each position, only one of the vectors contains data (i.e. the remaining are NA's).
#' 
#' Copied from \link[actel]{actel}
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

#' Remove a previously created transition layer
#' 
#' @export
#'
rmTransition <- function() {
  if (file.exists('rsp.transition.layer.RData')) {
    file.remove('rsp.transition.layer.RData')
    message("M: Transition layer removed")
  } else {
    message("M: No transition layer found in current working directory.")
  }
}
#' Identify and remove duplicated timestamps
#' 
#' @param input The detections to be checked
#' @param group The group being analysed, used solely for messaging purposes
#' 
#' @return The detections without duplicated timestamps
#' 
#' @keywords internal
#' 
checkDupTimestamps <- function(input, group, verbose = TRUE) {
  index <- which(duplicated(input$Timestamp))
  if (length(index) > 0) {
    input <- input[-index, ]
    if (verbose)
      warning(length(index), " individual detections were removed in group ", group," due to simultaneous detections at two receivers.", immediate. = TRUE, call. = FALSE)
  }
  return(input)
}

#' Performs a series of quality checks on the detection data.
#' 
#' @param input The detections list
#' @inheritParams dynBBMM
#' 
#' @return The detections which can be used for dbbmm
#' 
#' @keywords internal
#' 
checkGroupQuality <- function(input, verbose = TRUE) {
  if (attributes(input)$type == "group") {
    output <- lapply(names(input), function(i) {
      outp <- checkDupTimestamps(input = input[[i]], group = i, verbose = verbose)
      outp <- checkTrackPoints(input = outp, group = i, verbose = verbose)
      if (!is.null(outp)) {
        outp <- checkTrackTimes(input = outp, group = i, verbose = verbose)
        return(outp)
      } else {
        return(NULL)
      }
    })
    if (all(link <- unlist(lapply(output, is.null)))) {
      stop("All detection data failed to pass the quality checks for dBBMM implementation. Aborting.\n", call. = FALSE)
    }
    names(output) <- names(input)
    output <- output[!link]
    attributes(output)$type <- "group"
    return(output)
  }

  if (attributes(input)$type == "timeslot") {
    output <- lapply(names(input), function(g) {
      recipient <- lapply(names(input[[g]]), function(i) {
        aux <- checkDupTimestamps(input = input[[g]][[i]], group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
        aux <- checkTrackPoints(input = aux, group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
        if (!is.null(aux)) {
          aux <- checkTrackTimes(input = aux, group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
          return(aux)
        } else {
          return(NULL)
        }
      })
      names(recipient) <- names(input[[g]])
      return(recipient[!unlist(lapply(recipient, is.null))])
    })
    if (all(link <- unlist(lapply(output, length)) == 0))
      stop("All detection data failed to pass the quality checks for dBBMM implementation. Aborting.\n", call. = FALSE)
    names(output) <- names(input)
    output <- output[!link]
    attributes(output)$type <- "timeslot"
    return(output)
  }
}

#' Exclude tracks with less than 8 detections
#' 
#' @inheritParams checkDupTimestamps
#' 
#' @return The detections for tracks with more than 8 detections
#' 
#' @keywords internal
#' 
checkTrackPoints <- function(input, group, verbose = TRUE) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, nrow)) > 8
  if (all(!link)) {
    if (verbose)
      warning("ALL tracks in group ", group, " have less than eight detections. Removing group from analysis.", immediate. = TRUE, call. = FALSE)    
    return(NULL)
  } else {
    output <- tracks[link]
    if (verbose && length(tracks) > length(output))
      warning(length(tracks) - length(output), " track(s) in group ", group, " have less than eight detections and will not be used.", immediate. = TRUE, call. = FALSE)
    return(do.call(rbind.data.frame, output))
  }
}

#' Exclude tracks shorter than 30 minutes:
#' 
#' @inheritParams checkDupTimestamps
#' 
#' @return The detections for tracks longer than 30 minutes
#' 
#' @keywords internal
#' 
checkTrackTimes <- function(input, group, verbose = TRUE) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, function(x) {
    as.numeric(difftime(x$Timestamp[[nrow(x)]], x$Timestamp[[1]], units = "min"))
  })) >= 30
  if (all(!link)) {
    if (verbose)
      warning("ALL tracks in group ", group, " are shorter than 30 minutes. Removing group from analysis.", immediate. = TRUE, call. = FALSE)    
    return(NULL)
  } else {
    output <- tracks[link]
    if (verbose && length(tracks) > length(output))
      warning(sum(!link), " track(s) in group ", group, " are shorter than 30 minutes and will not be used.", immediate. = TRUE, call. = FALSE)
    return(do.call(rbind.data.frame, output))
  }
}


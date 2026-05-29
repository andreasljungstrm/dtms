#' Summarize transition probabilities
#'
#' @description
#' Provides several summary statistics on transition probabilities.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param fromvar Character (optional), name of variable with starting state in `probs`. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in `probs`. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in `probs`. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in `probs`. Default is `P`.
#' @param digits Numeric (optional), number of digits to return, default is 6.
#' @param format Character (optional), show results in decimal format or percentage, either `decimal` or `percent`. Default is `decimal`.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return A data frame
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' fit <- dtms_fit(data=estdata)
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' summary(probs)

dtms_probs_summary <- function(probs,
                               fromvar="from",
                               tovar="to",
                               timevar="time",
                               Pvar="P",
                               digits=4,
                               format="decimal",
                               sep="_") {

  # Get short state names
  probs[[fromvar]] <- dtms_getstate(probs[[fromvar]], sep=sep)
  probs[[tovar]]   <- dtms_getstate(probs[[tovar]],   sep=sep)

  # Aggregate in a single pass, then convert to data.frame for downstream ops
  dt     <- data.table::as.data.table(probs)
  result <- as.data.frame(dt[, .(
    MIN    = min(.SD[[1L]]),
    MAX    = max(.SD[[1L]]),
    MEDIAN = stats::median(.SD[[1L]]),
    MEAN   = mean(.SD[[1L]])
  ), by=c(fromvar, tovar), .SDcols=c(Pvar)])

  # Add time values for min and max
  result$MINtime <- probs[[timevar]][match(result$MIN, probs[[Pvar]])]
  result$MAXtime <- probs[[timevar]][match(result$MAX, probs[[Pvar]])]

  # Reorder columns to match original layout
  result <- result[, c(fromvar, tovar, "MIN", "MINtime", "MAX", "MAXtime", "MEDIAN", "MEAN")]

  # Order result
  ordering <- order(result[[fromvar]], result[[tovar]])
  result <- result[ordering,]

  # Rounding
  result[,c("MIN","MAX","MEDIAN","MEAN")] <-
    round(result[,c("MIN","MAX","MEDIAN","MEAN")],digits=digits)

  # For printing
  if(format=="percent") {
    result[,c("MIN","MAX","MEDIAN","MEAN")] <-
      result[,c("MIN","MAX","MEDIAN","MEAN")]*100
  }

  # Return
  return(result)

}

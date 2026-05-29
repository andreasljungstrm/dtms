#' Left censoring, right censoring, and gaps in data
#'
#' @description
#' This function provides an overview of censoring and gaps in the data. It can
#' do so in several ways: by providing counts of units with left censoring,
#' right censoring, and gaps; by providing a cross-tabulation of the number of
#' units with left censoring and/or right censoring and/or gaps; and by
#' returning a data frame with added indicators on censoring and gaps.
#'
#' @details
#' Added variables can be at the unit level or at the observation level. This
#' is controlled by the argument "addtype". If it is set to "id" then the unit
#' level is used. In this case the added variables are the same for
#' each observation  of a unit. For instance, if a unit experiences any gap,
#' then the added variable has the value TRUE for all observations of that unit.
#' If "addtype" is set to "obs" the observation level is used and the indicators
#' are only set to TRUE if they apply to a specific observation. For instance,
#' if a unit experience right censoring, only the last observation will have
#' TRUE as the value for the right-censoring indicator; i.e., showing that after
#' this last observation there is right censoring. This can be helpful for
#' analyses to understand censoring better.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param print Logical (optional), print counts? Default is TRUE.
#' @param printlong Logical (optional), print cross-tabulation? Default is FALSE.
#' @param add Logical (optional), add indicators to data set? Default is FALSE. If TRUE the data frame specified with \code{data} is returned with added columns.
#' @param addtype Character (optional), what type of information should be added if add=TRUE. Either "id" or "obs", see details. Default is "id".
#' @param varnames Character vector (optional), names of added variables if add=TRUE. Default is "c("LEFT","GAP","RIGHT")".
#'
#' @return Table or data frame.
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:19)
#' # Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Censoring
#' dtms_censoring(data=estdata,
#'                dtms=simple)

.dtms_censoring_unit_stats <- function(data, dtms, fromvar, tovar, timevar, idvar) {
  # Select only needed columns with fixed internal names (no get() overhead)
  dt <- data.table::as.data.table(data)[, c(idvar, timevar, tovar), with=FALSE]
  data.table::setnames(dt, c(".id", ".t", ".to"))

  tmin <- min(dtms$timescale)
  tmax <- max(dtms$timescale)
  step <- dtms$timestep
  absv <- paste(dtms$absorbing)

  data.table::setorderv(dt, c(".id", ".t"))

  # Vectorized boundary detection — no by= grouping needed.
  # After sorting by (id, t), group boundaries are where the id changes.
  next_id  <- data.table::shift(dt$.id, n=1L, type="lead")
  next_t   <- data.table::shift(dt$.t,  n=1L, type="lead")
  prev_id  <- data.table::shift(dt$.id, n=1L, type="lag")
  is_last  <- is.na(next_id) | (next_id != dt$.id)
  is_first <- is.na(prev_id) | (prev_id != dt$.id)
  # Time diff to next row within group (NA for last rows)
  td <- data.table::fifelse(is_last, NA_real_, next_t - dt$.t)

  # Per-row indicator columns (all fully vectorized)
  dt[, .left_f  := is_first & (dt$.t > tmin)]
  dt[, .gap_f   := !is_last & !(td %in% step)]   # gap to next obs != step
  dt[, .right_f := is_last & !(dt$.to %in% absv) & (dt$.t != tmax)]

  # One grouped aggregation using GForce-optimized sum (no .SD)
  stats <- dt[, .(
    left  = sum(.left_f)  > 0L,
    gap   = sum(.gap_f)   > 0L,
    right = sum(.right_f) > 0L
  ), by=.id]

  data.table::setnames(stats, ".id", idvar)
  stats
}

dtms_censoring <- function(data,
                           dtms,
                           fromvar="from",
                           tovar="to",
                           timevar="time",
                           idvar="id",
                           print=TRUE,
                           printlong=FALSE,
                           add=FALSE,
                           addtype="id",
                           varnames=c("LEFT","GAP","RIGHT")) {

  # Check dtms
  dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[[idvar]], data[[timevar]])
  data <- data[dataorder,]

  unit_stats <- .dtms_censoring_unit_stats(data, dtms, fromvar, tovar, timevar, idvar)
  left  <- unit_stats$left
  gap   <- unit_stats$gap
  right <- unit_stats$right

  # Print
  if(print) {
    cat("Units with left censoring: ",sum(left),"\n")
    cat("Units with gaps: ",sum(gap),"\n")
    cat("Units with right censoring: ",sum(right),"\n")
  }
  if(printlong) {
    cat("Cross tabulation:","\n")
    print(table(left,right,gap))
  }

  # Add indicators to data by ID
  if(add & addtype=="id") {
    tmp <- unit_stats
    names(tmp) <- c(idvar, varnames)
    data <- merge(data, tmp, by=idvar)
    return(data)
  }

  # Add indicators to data by observation
  if(add & addtype=="obs") {
    dt   <- data.table::as.data.table(data)
    tmin <- min(dtms$timescale)
    tmax <- max(dtms$timescale)
    step <- dtms$timestep
    absv <- paste(dtms$absorbing)
    # data is already sorted from the top of this function
    id_v  <- dt[[idvar]]
    t_v   <- dt[[timevar]]
    to_v  <- dt[[tovar]]
    next_id  <- data.table::shift(id_v, n=1L, type="lead")
    next_t   <- data.table::shift(t_v,  n=1L, type="lead")
    prev_id  <- data.table::shift(id_v, n=1L, type="lag")
    is_last  <- is.na(next_id) | (next_id != id_v)
    is_first <- is.na(prev_id) | (prev_id != id_v)
    td <- data.table::fifelse(is_last, NA_real_, next_t - t_v)
    dt[, (varnames[1]) := is_first & (t_v > tmin)]
    dt[, (varnames[2]) := !is_last & !(td %in% step)]
    dt[, (varnames[3]) := is_last & !(to_v %in% absv) & (t_v != tmax)]
    return(dt)
  }

}

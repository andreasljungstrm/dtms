#' Aggregate data
#'
#' @description
#' This function takes any data set and returns a new data frame which only
#' includes the unique rows from the original data set and indicates how
#' often these rows appear in the original data.
#'
#' @details
#' Currently, missing values are not supported and will be dropped; consider
#' using factors if you want to keep them. The variable provided with the
#' argument `idvar` will be dropped from the aggregated data. If `weights` is
#' specified, the counts will be placed in a variable with the same name. If
#' `countvar` is used, any existing variable in the original data with the
#' same name will be replaced.
#'
#' @param data Data frame.
#' @param weights Character (optional). Name of variable with weights.
#' @param idvar Character (optional). Name of variable in `data` with unit ID. Default is "id".
#' @param countvar Character (optional). Name of new variable in data which provides the counts. Default is "count".
#'
#' @returns An aggregated data frame
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' ## Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Aggregate
#' aggdata <- dtms_aggregate(estdata)
#' ## Fit model
#' fit <- dtms_fit(data=aggdata,
#'                 weights="count")

dtms_aggregate <- function(data,
                           weights=NULL,
                           idvar="id",
                           countvar="count") {

  # Transform to data frame, e.g., if tibble
  data <- data.table::as.data.table(data)

  # Get variable names without idvar (and weights if given)
  group_cols <- names(data)
  group_cols <- group_cols[group_cols != idvar]
  if(!is.null(weights)) group_cols <- group_cols[group_cols != weights]

  # Response column
  response_col <- if(is.null(weights)) {
    data.table::set(data, j=countvar, value=1L)
    countvar
  } else {
    weights
  }

  # Warning if missing values
  drops <- nrow(data) - nrow(stats::na.omit(data))
  if(drops > 0) warning(paste("Dropping", drops, "rows with missing values"))

  # Aggregate
  data <- stats::na.omit(data)
  tmp  <- data[, setNames(list(sum(get(response_col))), response_col), by=group_cols]

  # Return
  return(tmp)

}

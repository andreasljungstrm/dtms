#' Summarize data in transition format
#'
#' @description
#' Returns a data frame with number of observed transitions (column COUNT),
#' relative proportion of a transition relative to all transitions (column
#' PROP), and raw transition probabilities Pr(j|i) (column PROB).
#'
#' @param data Data frame, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable  in `data`with receiving state. Default is "to".
#' @param weights Character (optional), name of variable in `data` with weights. Default is NULL.
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
#' dtms_data_summary(estdata)

## Method
dtms_data_summary <- function(data,
                              dtms=NULL,
                              fromvar="from",
                              tovar="to",
                              weights=NULL) {

    # Convert to data.table
    data <- data.table::as.data.table(data)

    # Weights per transition
    if(is.null(weights)) data[, COUNT := 1L] else
      data.table::setnames(data, weights, "COUNT")

    # For handling of missing values
    data.table::set(data, i=which(is.na(data[[fromvar]])), j=fromvar, value="NA")
    data.table::set(data, i=which(is.na(data[[tovar]])),   j=tovar,   value="NA")

    # Aggregate
    result <- data[, .(COUNT=sum(COUNT)), by=c(fromvar, tovar)]

    # Order
    ordering <- order(result[[fromvar]], result[[tovar]])
    result <- result[ordering,]

    # Proportion
    N <- sum(result$COUNT)
    result[, PROP := COUNT / N]

    # Raw transition probabilities
    probs <- tapply(result$COUNT,
                    result[[fromvar]],
                    FUN=function(x) x/sum(x))
    probs <- unlist(probs)
    result[, PROB := probs]

    # If dtms is provided
    if(!is.null(dtms))  {

      # Get data frame
      newframe <- expand.grid(from=dtms$transient,to=c(dtms$transient,dtms$absorbing))

      # Result
      newresult <- merge(newframe,result,all=TRUE)

      # Replace missings with zero
      newresult[is.na(newresult$COUNT),c("COUNT","PROP","PROB")] <- rep(0,3)

      # Warning
      if(any(!result[[fromvar]]%in%c(dtms$transient,dtms$absorbing))) warning("Some fromvar values not in state space")
      if(any(!result[[tovar]]%in%c(dtms$transient,dtms$absorbing))) warning("Some tovar values not in state space")

    }

    # Return
    return(result)
}

#' Estimate (un)constrained discrete-time multistate model
#'
#' @description
#' This function estimates a (un)constrained discrete-time multistate model
#' using multinomial logistic regression.
#'
#' @details
#' The argument `data` takes a data set in transition format. The model formula
#' can either be specified by using the argument `formula`. Alternatively, it
#' can be specified with the arguments `fromvar`, `tovar`, and `controls`. These
#' are used if `formula` is not specified. `fromvar` takes the name of the
#' variable with the starting state as a character string, `tovar` the same for
#' the receiving state, and `controls` is an optional vector of control
#' variables. `fromvar` and `tovar` have default values which match other
#' functions of this package, making them a convenient alternative to `formula`
#' (see example).
#'
#' If `full=TRUE` a fully interacted model will be estimated in which each
#' control variable is interacted with all starting states. This is equivalent
#' to a full or unconstrained multistate model in which each transition is a
#' regression equation.
#'
#' The argument `package` is used choose the package used for estimation.
#' Currently, `VGAM` (default), `nnet`, and `mclogit` are supported.
#' The functions used for estimation are, respectively, `vgam`, `multinom`,
#' and `mblogit`. Arguments for these functions are passed via `...`.
#'
#' The argument `reference` sets the reference category for the multinomial
#' logistic regression. Weights for the regression can be passed via the
#' arguments `weights`. See the documentation of the package and function
#' used for estimation for details.
#'
#' When `cores > 1`, the data is split by starting state and one model is
#' fitted per starting state in parallel. This is equivalent to the
#' unconstrained model (`full=TRUE`) and returns an object of class
#' `dtms_multifit`. In this case, the `full` argument is ignored. Custom
#' `formula` arguments are not supported with `cores > 1`.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param full Logical (optional), estimate fully interacted model? Default is FALSE.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param cores Numeric (optional), number of cores for parallel per-state estimation. When `cores > 1`, one model is fitted per starting state in parallel (equivalent to unconstrained model). Default is 1.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Returns an object with class depending on the package used, or a
#'   `dtms_multifit` list of per-state models when `cores > 1`.
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)

dtms_fit <- function(data,
                     controls=NULL,
                     formula=NULL,
                     weights=NULL,
                     fromvar="from",
                     tovar="to",
                     reference=1,
                     package="VGAM",
                     full=FALSE,
                     cores=1,
                     ...) {

  # Parallel per-state fitting
  if(cores > 1) {

    if(!is.null(formula)) {
      warning("Custom formula not supported with cores > 1; fitting single model.")
      cores <- 1L
    } else {

      if(full) message("full=TRUE is ignored when cores > 1; per-state fitting produces the equivalent unconstrained model.")

      require(package, character.only=TRUE, quietly=TRUE)

      # Formula without fromvar since data will be split by starting state
      sub_formula <- dtms_formula(controls=controls,
                                  fromvar=NULL,
                                  tovar=tovar,
                                  full=FALSE)

      # Unique starting states present in the data
      from_states <- unique(data[, fromvar])

      # Extract weights vector before parallel loop
      if(!is.null(weights)) wts_all <- data[, weights] else wts_all <- NULL

      # Factor-encode the full tovar so all levels are known
      data[, tovar] <- as.factor(data[, tovar])
      all_to_levels <- levels(data[, tovar])

      `%dopar%` <- foreach::`%dopar%`
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)

      models <- foreach::foreach(state=from_states,
                                 .packages=c(package)) %dopar% {

        state_rows <- which(data[, fromvar] == state)
        sub_data <- data[state_rows, , drop=FALSE]
        sub_data[, tovar] <- factor(sub_data[, tovar], levels=all_to_levels)
        sub_data[, tovar] <- stats::relevel(sub_data[, tovar], ref=reference)
        environment(sub_formula) <- environment()

        if(!is.null(wts_all)) sub_wts <- wts_all[state_rows] else sub_wts <- NULL

        if(package == "VGAM") {
          VGAM::vgam(formula=sub_formula,
                     family=VGAM::multinomial(refLevel=reference),
                     data=sub_data,
                     weights=sub_wts,
                     ...)
        } else if(package == "nnet") {
          nnet::multinom(formula=sub_formula,
                         data=sub_data,
                         weights=sub_wts,
                         ...)
        } else if(package == "mclogit") {
          mclogit::mblogit(formula=sub_formula,
                           data=sub_data,
                           weights=sub_wts,
                           ...)
        }
      }

      parallel::stopCluster(cluster)

      names(models) <- from_states
      attr(models, "package")   <- package
      attr(models, "fromvar")   <- fromvar
      attr(models, "tovar")     <- tovar
      attr(models, "reference") <- reference
      attr(models, "controls")  <- controls
      class(models) <- c("dtms_multifit", "list")
      return(models)

    }
  }

  # Single-model path (original implementation)

  # Require package used for estimation (requireNamespace does not help here)
  require(package,character.only=TRUE,quietly=TRUE)

  # Build formula if not specified
  if(is.null(formula)) formula <- dtms_formula(controls=controls,
                                               fromvar=fromvar,
                                               tovar=tovar,
                                               full=full)

  # Make sure environment for formula is correct (ugh)
  environment(formula) <- environment()

  # Get weights if specified
  if(!is.null(weights)) weights <- data[,weights]

  # Factors (needed by most packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)

  # VGAM
  if(package=="VGAM") {

    # Estimate
    fitted <- VGAM::vgam(formula=formula,
                         family=VGAM::multinomial(refLevel=reference),
                         data=data,
                         weights=weights,
                         ...)
  }

  #nnet
  if(package=="nnet") {

    # Estimate
    fitted <- nnet::multinom(formula=formula,
                             data=data,
                             weights=weights,
                             ...)

  }

  #mclogit
  if(package=="mclogit") {

    # Estimate
    fitted <- mclogit::mblogit(formula=formula,
                               data=data,
                               weights=weights,
                               ...)

  }


  # Return results
  return(fitted)

}

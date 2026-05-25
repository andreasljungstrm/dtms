#' Predict transition probabilities
#'
#' @description
#' `dtms_transitions()` predicts transition probabilities based on a model
#' estimated with `dtms_fit()`.
#'
#' @details
#' Depending on the model specification, the prediction of transition
#' probabilities will require values for predictor variables which can be
#' specified with the argument `controls`. This is done using a named list
#' where each entry name must correspond to a variable name in the model.
#' For time-constant variables, each list entry is of length one and provides
#' a value for the corresponding time-constant variable. For time-varying
#' variables, each entry must have the length of the time scale minus one, and
#' provide a value for each (potential) transition in the model; i.e., starting
#' from time t=0, starting from time t=1, etc., until time t=T-1. Alternatively,
#' it can be of the same length as the time scale; in this case, the last value
#' is dismissed.
#'
#' If `vcov=TRUE` the full variance-covariance matrix of the transition
#' probabilities will be returned instead of the transition probabilities. If
#' `ci=TRUE`, confidence intervals will be returned. Note that the calculation
#' uses a normal approximation and results below 0 or above 1 are possible.
#'
#' The argument `dropvar` controls whether the covariate values used for
#' prediction are dropped. If `FALSE` each row of the resulting data frame will
#' have the covariate values which were used to predict the corresponding
#' probability.
#'
#' @param model Model estimated with \code{dtms_fit}, or a \code{dtms_multifit}
#'   object from parallel per-state fitting.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param controls List (optional) with values for predictors (see details).
#' @param se Logical (optional), return standard errors of predicted probabilities. Default is `TRUE`.
#' @param vcov Logical (optional), return variance-covariance matrix of predicted probabilities. Default is `FALSE`.
#' @param ci Logical (optional), return confidence intervals? See details. Default is FALSE.
#' @param alpha Numeric (optional), if ci=TRUE, what confidence level is used? Default is 0.05.
#' @param dropvar Logical (optional), should covariate values used for prediction be returned (see details). Default is `TRUE`.
#' @param fromvar Character (optional), name of variable with starting state in the returned data frame. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in the returned data frame. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in the returned data frame. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in the returned data frame. Default is `P`.
#'
#' @return A data frame with transition probabilities.
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
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)

dtms_transitions <- function(model,
                             dtms,
                             controls=NULL,
                             dropvar=TRUE,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P",
                             se=TRUE,
                             vcov=FALSE,
                             ci=FALSE,
                             alpha=0.05) {

  # Check
  dtms_proper(dtms)

  # Adjust time scale (transitions in the model)
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]
  ntime <- length(timescale_reduced)

  # Get full state space
  all_states <- paste(c(dtms$transient,dtms$absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=dtms$transient,
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)

  # Get names right
  names(model_frame) <- c(fromvar,timevar)

  # Deal with controls
  varnames <- names(controls)
  for(var in varnames) {

    # Get values
    value <- controls[[var]]
    nvalue <- length(value)

    # Check
    if(!nvalue%in%c(1,ntime,ntime+1)) stop("Wrong number of time-varying values")

    # Act depending on how many values
    if(nvalue==1) model_frame[var] <- value

    if(nvalue==ntime) {
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }

    if(nvalue==ntime+1) {
      value <- value[-nvalue]
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }


  }

  # Predict for dtms_multifit (parallel per-state models)
  if(inherits(model, "dtms_multifit")) {

    pkg <- attr(model, "package")

    # Initialise all state columns to 0
    model_frame[, all_states] <- 0

    for(state in names(model)) {

      sub_model  <- model[[state]]
      state_rows <- which(model_frame[, fromvar] == state)

      if(length(state_rows) == 0) next

      sub_frame <- model_frame[state_rows, , drop=FALSE]

      if(inherits(sub_model, c("vgam", "mclogit"))) {
        sub_pred <- stats::predict(sub_model, sub_frame, "response")
      } else if(inherits(sub_model, "nnet")) {
        sub_pred <- stats::predict(sub_model, sub_frame, "probs")
      }

      seen_states <- colnames(sub_pred)
      model_frame[state_rows, seen_states] <- sub_pred[, seen_states, drop=FALSE]
    }

  }

  # Predict for vgam/mclogit
  if(inherits(model,c("vgam","mclogit"))) {

    # Actual prediction
    predicted <- stats::predict(model,model_frame,"response")
    predicted_states <- colnames(predicted)

    # Check if all states are there, if not place 0
    if(!all(all_states%in%predicted_states)) {
      warning("Some states defined via dtms do not appear in data, might cause issues with SEs")
      which_states <- all_states[which(!all_states%in%predicted_states)]
      n_states <- length(which_states)
      tmp <- matrix(data=0,ncol=n_states,nrow=dim(predicted)[1])
      colnames(tmp) <- which_states
      predicted <- cbind(predicted,tmp)
      model_frame[,all_states] <- predicted[,all_states]
    } else model_frame[,all_states] <- predicted[,all_states]
  }

  # Predict for nnet
  if(inherits(model,"nnet")) {
    predicted <- stats::predict(model,model_frame,"probs")
    predicted_states <- colnames(predicted)
    if(!all(all_states%in%predicted_states)) {
      warning("Some states defined via dtms do not appear in data, not supported with nnet")
    }
    model_frame[,all_states] <- predicted[,all_states]
  }

  # Values of starting state
  model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

  # Reshape
  model_frame <- stats::reshape(model_frame,
                                varying=all_states,
                                idvar=fromvar,
                                timevar=tovar,
                                times=all_states,
                                direction="long",
                                v.names=Pvar)

  # SE/CI/vcov
  if(se|vcov|ci) {

    # Simplify starting state (needed for model.matrix below)
    model_frame[,fromvar] <- dtms_simplify(model_frame)$from

    # --- dtms_multifit: per-state delta-method SE ---
    if(inherits(model, "dtms_multifit")) {

      pkg <- attr(model, "package")
      nstates <- length(all_states)
      nprobs_total <- nrow(model_frame)

      se_vec <- numeric(nprobs_total)

      if(vcov) {
        Vp_full <- matrix(0, nrow=nprobs_total, ncol=nprobs_total)
      }

      for(state in names(model)) {

        sub_model  <- model[[state]]
        state_rows <- which(model_frame[, fromvar] == state)

        if(length(state_rows) == 0) next

        sub_frame  <- model_frame[state_rows, , drop=FALSE]
        nprobs_sub <- length(state_rows)

        # Coefficients
        if(inherits(sub_model, "mclogit")) {
          C_sub       <- stats::coef(sub_model)
          Cstates_sub <- dtms_getstate(names(C_sub), sep="~")
          Cstates_sub <- unique(Cstates_sub)
          C_sub       <- matrix(data=C_sub,
                                ncol=length(Cstates_sub),
                                byrow=TRUE)
        }

        if(inherits(sub_model, "vgam")) {
          C_sub       <- stats::coef(sub_model)
          Cstates_sub <- sub_model@extra$colnames.y[-sub_model@extra$use.refLevel]
          C_sub       <- matrix(C_sub,
                                ncol=length(Cstates_sub),
                                byrow=TRUE)
        }

        if(inherits(sub_model, "nnet")) {
          C_sub       <- stats::coef(sub_model)
          Cstates_sub <- rownames(C_sub)
          C_sub       <- t(C_sub)
        }

        # vcov of sub-model coefficients
        Vml_sub <- stats::vcov(sub_model)

        if(inherits(sub_model, "mclogit")) {
          ordernames_sub <- sort(rownames(Vml_sub))
          Vml_sub        <- Vml_sub[ordernames_sub, ordernames_sub]
          Vnames_sub     <- dtms_getstate(rownames(Vml_sub), sep="~")
        }

        if(inherits(sub_model, "vgam")) {
          nicenames_sub <- dtms_getstate(rownames(Vml_sub), sep=":")
          nicenames_sub <- unique(nicenames_sub)
          nicenames_sub <- sort(dtms_combine(Cstates_sub, nicenames_sub, sep=":"))

          newnames_sub <- unlist(lapply(strsplit(colnames(Vml_sub), split=":"),
                                        function(x) paste0(x[2],":",x[1])))
          colnames(Vml_sub) <- rownames(Vml_sub) <- newnames_sub
          ordernames_sub    <- sort(rownames(Vml_sub))
          Vml_sub           <- Vml_sub[ordernames_sub, ordernames_sub]
          colnames(Vml_sub) <- rownames(Vml_sub) <- nicenames_sub

          Vnames_sub <- dtms_getstate(rownames(Vml_sub), sep=":")
        }

        if(inherits(sub_model, "nnet")) {
          Vnames_sub <- dtms_getstate(rownames(Vml_sub), sep=":")
        }

        ncoef_sub <- dim(Vml_sub)[1]

        # Model matrix for this sub-frame
        form_sub <- stats::formula(sub_model)
        mm_sub   <- stats::model.matrix(object=form_sub, data=sub_frame)

        # Scores (denominator for predicted prob)
        dscores_sub            <- matrix(data=1, ncol=nstates, nrow=nprobs_sub)
        colnames(dscores_sub)  <- sort(all_states)
        dscores_sub[, Cstates_sub] <- exp(mm_sub %*% C_sub)

        # Full score
        fullscore_sub <- rowSums(dscores_sub)

        Z_sub <- matrix(data=fullscore_sub, nrow=nprobs_sub, ncol=ncoef_sub)

        N_sub <- stats::model.matrix(object=~to, data=sub_frame)
        N_sub[,1] <- N_sub[,1] - rowSums(N_sub[,-1])
        N_sub <- rowSums(N_sub * dscores_sub)
        N_sub <- matrix(data=N_sub, nrow=nprobs_sub, ncol=ncoef_sub)

        varvalues_sub <- do.call("cbind", replicate(length(Cstates_sub), mm_sub, simplify=FALSE))
        Zdash_sub     <- varvalues_sub * dscores_sub[, Vnames_sub]

        Ndash_sub <- outer(sub_frame[, tovar], Vnames_sub, FUN=`==`)
        Ndash_sub <- Ndash_sub * N_sub * varvalues_sub

        G_sub <- (Ndash_sub * Z_sub - N_sub * Zdash_sub) / (Z_sub^2)

        if(inherits(sub_model, "mclogit")) colnames(G_sub) <- dtms_combine(Cstates_sub, colnames(mm_sub), sep="~")
        if(inherits(sub_model, "vgam"))    colnames(G_sub) <- dtms_combine(Cstates_sub, colnames(mm_sub), sep=":")
        if(inherits(sub_model, "nnet"))    colnames(G_sub) <- dtms_combine(Cstates_sub, colnames(mm_sub), sep=":")

        G_sub <- G_sub[, rownames(Vml_sub)]

        Vp_sub <- G_sub %*% Vml_sub %*% t(G_sub)

        se_vec[state_rows] <- sqrt(diag(Vp_sub))

        if(vcov) Vp_full[state_rows, state_rows] <- Vp_sub
      }

      # Return vcov matrix
      if(vcov) return(Vp_full)

      # Restore full starting state
      model_frame[, fromvar] <- paste(model_frame[, fromvar],
                                      model_frame[, timevar],
                                      sep=dtms$sep)

      if(se) model_frame$se <- se_vec

      if(ci) {
        z     <- stats::qnorm(1 - alpha/2)
        model_frame$cilow <- model_frame[, Pvar] - z * se_vec
        model_frame$ciup  <- model_frame[, Pvar] + z * se_vec
      }

    } else {

      # --- Single model: original delta-method SE ---

      # Coefficients
      if(inherits(model,"mclogit")) {
        C <- stats::coef(model)
        Cstates <- dtms_getstate(names(C),sep="~")
        Cstates <- unique(Cstates)
        C <- matrix(data=C,
                    ncol=length(all_states)-1,
                    byrow=TRUE)
      }

      if(inherits(model,"vgam")) {
        C <- stats::coef(model)
        C <- matrix(C,
                    ncol=length(all_states)-1,
                    byrow=TRUE)
        Cstates <- model@extra$colnames.y[-model@extra$use.refLevel]
      }

      if(inherits(model,"nnet")) {
        C <- stats::coef(model)
        Cstates <- rownames(C)
        C <- t(C)
      }

      # vcov of coefficients
      Vml <- stats::vcov(model)

      if(inherits(model,"mclogit")) {
        ordernames <- sort(rownames(Vml))
        Vml <- Vml[ordernames,ordernames]
        Vnames <- dtms_getstate(rownames(Vml),sep=c("~"))
      }

      if(inherits(model,"vgam")) {

        # Get nice names (assigned below)
        nicenames <- dtms_getstate(rownames(Vml),sep=c(":"))
        nicenames <- unique(nicenames)
        nicenames <- sort(dtms_combine(Cstates,nicenames,sep=":"))

        # Reorder
        newnames <- unlist(lapply(strsplit(colnames(Vml),split=":"),function(x) paste0(x[2],":",x[1])))
        colnames(Vml) <- rownames(Vml) <- newnames
        ordernames <- sort(rownames(Vml))
        Vml <- Vml[ordernames,ordernames]

        # Assign nice names
        colnames(Vml) <- rownames(Vml) <- nicenames

        # State names
        Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))
      }

      if(inherits(model,"nnet")) Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))

      # Number of probabilities, coefs, states
      nprobs <- dim(model_frame)[1]
      ncoef <- dim(Vml)[1]
      nstates <- length(all_states)

      # Model matrix
      form <- stats::formula(model)
      mm <- stats::model.matrix(object=form,data=model_frame)

      # Scores (denominator for predicted prob)
      dscores <- matrix(data=1,
                        ncol=nstates,
                        nrow=nprobs)
      colnames(dscores) <- sort(all_states)
      dscores[,Cstates] <- exp(mm%*%C)

      # Full score (numerator for predicted prob)
      fullscore <- rowSums(dscores)

      # Parts of full derivative (n'*z-n*z')/z^2
      Z <- matrix(data=fullscore,
                  nrow=nprobs,
                  ncol=ncoef)

      N <- stats::model.matrix(object=~to,data=model_frame)
      N[,1] <- N[,1]-rowSums(N[,-1])
      N <- rowSums(N*dscores)
      N <- matrix(data=N,
                  nrow=nprobs,
                  ncol=ncoef)

      varvalues <- do.call("cbind",replicate(nstates-1,mm,simplify=FALSE))
      Zdash <- varvalues*dscores[,Vnames]

      Ndash <- outer(model_frame$to,Vnames,FUN=`==`)
      Ndash <- Ndash*N*varvalues

      # Matrix of derivatives
      G <- (Ndash*Z-N*Zdash)/(Z^2)

      # Re-order
      if(inherits(model,"mclogit")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep="~")
      if(inherits(model,"vgam")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")
      if(inherits(model,"nnet")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")

      G <- G[,rownames(Vml)]

      # Vcov matrix
      Vp <- G%*%Vml%*%t(G)

      # Return vcov matrix
      if(vcov) return(Vp)

      # Full starting state
      model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

      # SE
      if(se) model_frame$se <- sqrt(diag(Vp))

      # CI
      if(ci) {
        z <- (1-alpha/2)
        z <- stats::qnorm(z)
        error <- sqrt(diag(Vp))
        model_frame$cilow <- model_frame[,Pvar]-z*error
        model_frame$ciup <- model_frame[,Pvar]+z*error
      }

    }

  }

  # Values of receiving state (state name + time)
  rightrows <- model_frame[,tovar]%in%dtms$transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+dtms$timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=dtms$sep)

  # Drop row names
  rownames(model_frame) <- NULL

  # Drop covariate values for prediction
  if(dropvar) {
    model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se","cilow","ciup")]
  }

  # Class
  class(model_frame) <- c("dtms_probs","data.frame")

  # Return
  return(model_frame)

}

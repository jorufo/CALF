#' @import data.table

calculateCalf <- function(resp, pred, maxPred, score, margin, verbose) {
  
  ## Parameter checks
  chkNum <- function(x) is.numeric(type.convert(x))
  rl <- length(resp)
  stopifnot(class(score) == "character", length(score) == 1)
  if (!is.null(margin)) {
    stopifnot(class(margin) == "numeric", length(margin) == 1)
  }
  stopifnot(class(verbose) == "logical", length(verbose) == 1)
  stopifnot(class(maxPred) == "numeric", length(maxPred) == 1)
  stopifnot(class(pred) %in% c("matrix", "data.frame", "data.table"))
  stopifnot(class(resp) %in% c("logical", "numeric"))
  stopifnot(rl == nrow(pred))
  chkNum <- function(x) is.numeric(type.convert(x))
  if (!all(apply(pred, 2, chkNum))) {
    stop("Some of the predictors in 'pred' are not coercible to numeric.")
  }
  if (!is.logical(resp) & score != "mse") {
    stop("'resp' assumed continuous; 'score' must be 'mse' -- see ?calf")
  }
  
  sFunc <- switch(score,
                  "auc" = calfScoreAUC,
                  "pcc" = calfScorePCC,
                  "wts" = calfScoreWTS,
                  FALSE)
  
  if (!is.function(sFunc)) stop("Invalid 'score' parameter -- see ?calf.")
  
  if (verbose) {
    stime <- Sys.time()
    msg <- sprintf("%17s  %9s  %8s  %5s  %6s",
                   "Elapsed (d:h:m:s)", "Iteration", 
                   "Score", "Param", "Weight")
    cat(msg, "\n")
  }
  
  if (!is.data.table(pred)) pred <- as.data.table(pred)
  
  s <- 0
  v <- vector("numeric", rl)
  i <- 1L
  np <- ncol(pred)
  pNames <- names(pred)
  
  nPred <- min(maxPred, np)
  pVec <- rep(NA_integer_, length = nPred) ## store the parameter index
  wVec <- rep(NA_integer_, length = nPred) ## store the weight
  tVec <- rep(NA_real_,    length = nPred) ## store the score @step "training"
  
  repeat {
    
    sVec <- c(sapply(v + pred, sFunc, resp = resp), 
              sapply(v - pred, sFunc, resp = resp))
    
    newPred <- which.max(sVec)
    
    ## Check if the best improvement in score surpasses the user-defined margin
    if (!is.null(margin) && !sVec[newPred] - s > margin) break
    
    neg <- newPred > np
    newPred <- newPred - np*neg
    
    if (i == 1L) {
      w <- sign(cor(pred[[newPred]], resp)) 
    } else {
      w <- 1L - 2L*neg
    }
    
    pVec[i] <- match(names(newPred), pNames)
    wVec[i] <- w
    tVec[i] <- sVec[[newPred + np*neg]]
    
    v <- v + w*pred[[newPred]]
    pred[[newPred]] <- NULL
    np <- ncol(pred)
    s <- tVec[i]
    
    if (verbose) {
      etime <- Sys.time() - stime
      units(etime) <- "secs"
      etime <- unclass(etime)
      etime <- sprintf("%02d:%02d:%02d:%02d", 
                       etime %/% 86400,  
                       etime %% 86400 %/% 3600,  
                       etime %% 3600 %/% 60,
                       etime %% 60 %/% 1)
      msg <- sprintf("%17s  %9s  %8s  %5s  %6s", 
                     etime, i, signif(tVec[i], 5), pVec[i], wVec[i])
      cat(msg, "\n")
    }
    
    i <- i + 1L
    if (i > nPred) break
    
  }
  
  res <- list(predInd  = pVec,
              predName = pNames[pVec],
              weight   = wVec,
              trainVec = tVec,
              score    = score,
              maxPred  = maxPred,
              margin   = margin)
  
  class(res) <- "calf"
  
  res
  
}

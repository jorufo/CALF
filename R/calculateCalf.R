#' @import data.table

calculateCalf <- function(x, y, maxX, score, margin, verbose) {
  
  ## Parameter checks
  chkNum <- function(x) is.numeric(type.convert(x))
  yl <- length(y)
  stopifnot(class(score) == "character", length(score) == 1)
  if (!is.null(margin)) {
    stopifnot(class(margin) == "numeric", length(margin) == 1)
  }
  stopifnot(class(verbose) == "logical", length(verbose) == 1)
  stopifnot(class(maxX) == "numeric", length(maxX) == 1)
  stopifnot(class(x) %in% c("matrix", "data.frame", "data.table"))
  stopifnot(class(y) %in% c("logical", "numeric"))
  stopifnot(yl == nrow(x))
  chkNum <- function(x) is.numeric(type.convert(x))
  if (!all(apply(x, 2, chkNum))) {
    stop("Some of the predictors in 'x' are not coercible to numeric.")
  }
  if (!is.logical(y) & score != "mse") {
    stop("'y' assumed continuous; 'score' must be 'mse' -- see ?calf")
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
  
  hasNames <- !is.null(colnames(x))
  if (!is.data.table(x)) x <- as.data.table(x)
  
  s <- 0
  v <- vector("numeric", yl)
  i <- 1L
  nX <- ncol(x)
  xNames <- names(x)
  
  maxIter <- min(maxX, nX)
  xVec <- rep(NA_integer_, length = maxIter) ## store the parameter index
  wVec <- rep(NA_integer_, length = maxIter) ## store the weight
  tVec <- rep(NA_real_,    length = maxIter) ## store the score @step "training"
  
  repeat {
    
    sVec <- c(sapply(v + x, sFunc, y = y), 
              sapply(v - x, sFunc, y = y))
    
    newX <- which.max(sVec)
    
    ## Check if the best improvement in score surpasses the user-defined margin
    if (!is.null(margin) && !sVec[newX] - s > margin) break
    
    neg <- newX > nX
    newX <- newX - nX*neg
    
    if (i == 1L) {
      w <- sign(cor(x[[newX]], y)) 
    } else {
      w <- 1L - 2L*neg
    }
    
    xVec[i] <- match(names(newX), xNames)
    wVec[i] <- w
    tVec[i] <- sVec[[newX + nX*neg]]
    
    v <- v + w*x[[newX]]
    x[[newX]] <- NULL
    nX <- ncol(x)
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
                     etime, i, signif(tVec[i], 5), xVec[i], wVec[i])
      cat(msg, "\n")
    }
    
    i <- i + 1L
    if (i > maxIter) break
    
  }
  
  res <- list(xInd     = xVec,
              xName    = xNames[xVec],
              weight   = wVec,
              trainVec = tVec,
              score    = score,
              maxX     = maxX,
              margin   = margin)
  
  class(res) <- "calf"
  
  res
  
}

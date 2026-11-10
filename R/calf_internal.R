calf_score <- function(target, values, targetVector, optimize) {
  if (targetVector == "nonbinary") {
    value <- suppressWarnings(cor(target, values))
    return(if (is.finite(value) && value > 0) 1 / value else Inf)
  }
  caseVar <- values[target == 1]
  ctrlVar <- values[target == 0]
  if (optimize == "auc") {
    value <- compute.auc(caseVar, ctrlVar)
    return(if (is.finite(value) && value > 0) 1 / value else Inf)
  }
  value <- tryCatch(suppressWarnings(t.test(caseVar, ctrlVar)$p.value),
                    error = function(e) NA_real_)
  if (is.finite(value)) value else Inf
}

calf_internal <- function(data,
                          nMarkers,
                          randomize  = FALSE,
                          proportion = NULL,
                          times,
                          targetVector = "binary",
                          optimize = "pval",
                          verbose = FALSE){
  
  
  # getting rid of global variable warning -------------------------- #
  x = NULL
  y = NULL
  refx = NULL
  refy = NULL
  
  if (targetVector == "nonbinary")
      optimize <- NULL
  
  # setting up some initial values -----------------------------------#
  if ((!is.data.frame(data) && !is.matrix(data)) || ncol(data) < 2L ||
      nrow(data) < 3L || any(!vapply(as.data.frame(data), is.numeric, logical(1)))) {
    stop("CALF ERROR: Data are not numeric. Please check that data were read in correctly.")
  }
  if (anyNA(data) || any(!is.finite(as.matrix(data))))
    stop("CALF ERROR: Data must contain only finite, non-missing values.")
  if (length(nMarkers) != 1L || !is.finite(nMarkers) ||
      nMarkers != as.integer(nMarkers) || nMarkers < 1L)
    stop("CALF ERROR: nMarkers must be a positive integer.")
  if (length(targetVector) != 1L || !targetVector %in% c("binary", "nonbinary"))
    stop("CALF ERROR: targetVector must be 'binary' or 'nonbinary'.")
  if (targetVector == "binary" && !all(data[, 1] %in% c(0, 1)))
    stop("CALF ERROR: Binary targets must contain only zero and one.")
  if (targetVector == "binary" && !optimize %in% c("pval", "auc"))
    stop("CALF ERROR: Binary optimize must be 'pval' or 'auc'.")
  
  
  nVars <- ncol(data) - 1
  dNeg  <- data[ ,2:ncol(data), drop = FALSE]
  dNeg  <- dNeg * - 1
  data  <- data.frame(data, dNeg, check.names = FALSE)
  
  if (nMarkers > nVars){
    stop(paste0("CALF ERROR: Requested number of markers is larger than the number of markers in data set. ",
                "Please revise this value or make sure your data were read in properly."))
  }
  
  if (randomize == TRUE) data[ ,1] <- sample(data[ ,1])
  
  if (!is.null(proportion)) {
    if (length(proportion) < 1L || length(proportion) > 2L ||
        any(!is.finite(proportion)) || any(proportion <= 0 | proportion > 1) ||
        (targetVector == "nonbinary" && length(proportion) != 1L))
      stop("CALF ERROR: Invalid sampling proportion.")
    
    if (targetVector == "binary") {
      
      ctrlRows  <- which(data[ ,1] == 0)
      caseRows  <- which(data[ ,1] == 1)
      
      # calculate number of case and control to keep
      
      if(length(proportion) == 2) {
        
        nCtrlKeep <- round(length(ctrlRows)*proportion[1], digits = 0)
        nCaseKeep <- round(length(caseRows)*proportion[2], digits = 0)
        
      } else {
        
        nCtrlKeep <- round(length(ctrlRows)*proportion, digits = 0)
        nCaseKeep <- round(length(caseRows)*proportion, digits = 0)
        
      }
      
      # sample randomly rows of case and control to keep, record rows to keep
      if (min(nCtrlKeep, nCaseKeep) < 2L)
        stop("CALF ERROR: At least two cases and controls must be retained.")
      keepRows <- c(sample(ctrlRows, nCtrlKeep), sample(caseRows, nCaseKeep))
      
      # subset original data to keep these rows
      data      <- data[keepRows, , drop = FALSE]
      
    } else {
      
      nDataKeep <- round(nrow(data)*proportion, digits = 0)
      if (nDataKeep < 3L)
        stop("CALF ERROR: At least three observations must be retained.")
      keepRows <- sample(seq_len(nrow(data)), nDataKeep)
      data <- data[keepRows, , drop = FALSE]
      
    }
    
  }
  
  real  <- data[ ,1]
  realMarkers <- data[ , 2:ncol(data), drop = FALSE]
  ctrl  <- data[data[ ,1] == 0, 2:ncol(data), drop = FALSE]
  case  <- data[data[ ,1] == 1, 2:ncol(data), drop = FALSE]
  if (targetVector == "binary" && (nrow(ctrl) < 2L || nrow(case) < 2L))
    stop("CALF ERROR: At least two cases and two controls are required.")
  indexNegPos <- rep(0, (nVars*2))
  # end of setting up some initial values ----------------------------#
  
  # initial loop to establish first optimal marker -------------------#
  allCrit <- numeric()
  for (i in 1:(nVars*2)){
    allCrit[i] <- calf_score(real, realMarkers[, i], targetVector, optimize)
  }
  if (!any(is.finite(allCrit)))
    stop("CALF ERROR: No marker has a valid optimization score.")
  
  # end of initial loop ----------------------------------------------#
  
  keepMarkers  <- names(realMarkers)[which.min(allCrit)]
  bestCrit     <- min(allCrit)
  keepIndex    <- which.min(allCrit)
  
  if (verbose == TRUE) {
    
    if(targetVector == "binary") {
      
      if (optimize == "pval"){
        cat("Selected:", keepMarkers,
            paste0("p value = ", round(bestCrit, digits = 15), "\n"))
      } else if (optimize == "auc"){
        cat("Selected:", keepMarkers,
            paste0("AUC = ", round((1/bestCrit), digits = 15), "\n"))
      }
      
    } else if (targetVector == "nonbinary") {
        cat("Selected:", keepMarkers,
            paste0("Correlation = ", round((1/bestCrit), digits = 15), "\n"))
    }
  }
  
  # second loop to add another marker --------------------------------#
  if (nMarkers != 1){
    allCrit  <- numeric()
    realPrev <- realMarkers[ ,keepIndex]
    casePrev <- case[ ,keepIndex]
    ctrlPrev <- ctrl[ ,keepIndex]
    for (i in 1:(nVars*2)){
      # Check the indicies and complement for the postive values, else the negative ones
      if( i >= 1 && i <= nVars ) {
        # ensure an index or its complement is not being used if already chosen 
        if (i != keepIndex && (nVars+i) != keepIndex){
          caseVar <- casePrev + case[ ,i]
          ctrlVar <- ctrlPrev + ctrl[ ,i]
          realVar <- realPrev + realMarkers[ ,i]
          crit <- calf_score(real, realVar, targetVector, optimize)
        } else {
          crit <- NA
        }
      } else if( i >= (nVars+1) && i <= 2*nVars) {
        # ensure an index or its complement is not being used if already chosen 
        if (i != keepIndex && (i-nVars) != keepIndex){
          caseVar <- casePrev + case[ ,i]
          ctrlVar <- ctrlPrev + ctrl[ ,i]
          realVar <- realPrev + realMarkers[ ,i]
          crit <- calf_score(real, realVar, targetVector, optimize)
        } else {
          crit <- NA
        }
      } else {
        crit <- NA  # Should really never get here
      }
      allCrit[i] <- crit
    }
    # end of second loop ----------------------------------------------#
    
    # An absent valid candidate is a normal termination condition.
    continue <- any(is.finite(allCrit)) &&
      bestCrit[length(bestCrit)] > min(allCrit, na.rm = TRUE)
    
    
    if (continue == TRUE){
      keepMarkers  <- append(keepMarkers, names(realMarkers)[which.min(allCrit)])
      bestCrit     <- append(bestCrit, min(allCrit, na.rm = TRUE))
      keepIndex    <- append(keepIndex, which.min(allCrit))
      
      if (length(keepMarkers) == nMarkers) continue <- FALSE
    }
    
    if (verbose == TRUE) {
      
      if(targetVector == "binary") {

        if (optimize == "pval"){
          cat("Selected:", keepMarkers[length(keepMarkers)],
              paste0("p value = ", round(bestCrit[length(bestCrit)], digits = 15), "\n"))
        } else if (optimize == "auc"){
          cat("Selected:", keepMarkers[length(keepMarkers)],
              paste0("AUC = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
        }
        
      } else if (targetVector == "nonbinary") {
        cat("Selected:", keepMarkers[length(keepMarkers)],
            paste0("Correlation = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
      }
    
    }
    
    # loop for third through nMarker ----------------------------------#
    while (continue == TRUE){
      allCrit  <- numeric()
      casePrev <- rowSums(case[ ,keepIndex, drop = FALSE])
      ctrlPrev <- rowSums(ctrl[ ,keepIndex, drop = FALSE])
      realPrev <- rowSums(realMarkers[ ,keepIndex, drop = FALSE])
      for (i in 1:(nVars*2)){
        if( i >= 1 && i <= nVars ) {
          if ( !(i %in% keepIndex) && !((nVars+i) %in% keepIndex) ){
            caseVar <- casePrev + case[ ,i]
            ctrlVar <- ctrlPrev + ctrl[ ,i]
            realVar <- realPrev + realMarkers[ ,i]
              crit <- calf_score(real, realVar, targetVector, optimize)
          } else {
            crit <- NA
          }
        } else if( i >= (nVars+1) && i <= 2*nVars) {
          if ( !(i %in% keepIndex) && !((i-nVars) %in% keepIndex) ){
            caseVar <- casePrev + case[ ,i]
            ctrlVar <- ctrlPrev + ctrl[ ,i]
            realVar <- realPrev + realMarkers[ ,i]
              crit <- calf_score(real, realVar, targetVector, optimize)
          } else {
            crit <- NA
          }
        } else {
          crit <- NA  # Should really never get here
        }
        
        allCrit[i] <- crit
      }
      continue <- any(is.finite(allCrit)) &&
        bestCrit[length(bestCrit)] > min(allCrit, na.rm = TRUE)

      if (continue == TRUE){
        keepMarkers  <- append(keepMarkers, names(realMarkers)[which.min(allCrit)])
        bestCrit     <- append(bestCrit, min(allCrit, na.rm = TRUE))
        keepIndex    <- append(keepIndex, which.min(allCrit))
        continue     <- bestCrit[length(bestCrit)] < bestCrit[length(bestCrit)-1]
        if (verbose == TRUE) {
          
          if(targetVector == "binary") {
            if (optimize == "pval"){
              cat("Selected:", keepMarkers[length(keepMarkers)],
                  paste0("p value = ", round(bestCrit[length(bestCrit)], digits = 15), "\n"))
            } else if (optimize == "auc"){
              cat("Selected:", keepMarkers[length(keepMarkers)],
                  paste0("AUC = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
            }
          } else if (targetVector == "nonbinary") {
            cat("Selected:", keepMarkers[length(keepMarkers)],
                paste0("Correlation = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
          }
          
        }
      }
      
      if (length(keepMarkers) == nMarkers) continue <- FALSE
    }
  }
  
  if (verbose == TRUE) cat("\n")
  
  indexNegPos[keepIndex] <- ifelse(keepIndex > nVars, -1, 1)
  finalIndex   <- ifelse(keepIndex <= nVars, keepIndex, keepIndex - nVars)
  finalMarkers <- data.frame(names(case)[finalIndex], indexNegPos[keepIndex], check.names = FALSE)
  names(finalMarkers) <- c("Marker","Weight")
  
  if (targetVector == "nonbinary" || optimize == "auc") {
    finalBestCrit <- 1 / bestCrit[length(bestCrit)]
  } else {
    finalBestCrit <- bestCrit[length(bestCrit)]
  }
  ## AUC -------------------------------------------------------------#
  # create function value for each individual
  if (targetVector == "binary"){
    scores <- round(rowSums(realMarkers[, keepIndex, drop = FALSE]), 8)
    labels <- data[, 1]
    auc <- compute.auc(scores[labels == 1], scores[labels == 0])
    # CALF reports the better of both score orientations.
    ranks <- rank(round(scores, 8), ties.method = "average")
    ncase <- sum(labels == 1)
    nctrl <- sum(labels == 0)
    caseAuc <- (sum(ranks[labels == 1]) - ncase * (ncase + 1) / 2) /
      (ncase * nctrl)
    if (caseAuc < 0.5) scores <- -scores
    ord <- order(scores, decreasing = TRUE)
    orderedScores <- scores[ord]
    orderedLabels <- labels[ord]
    # Tied observations enter at the same threshold.
    lastOfTie <- !duplicated(orderedScores, fromLast = TRUE)
    rocData <- data.frame(
      x = c(0, cumsum(orderedLabels == 0)[lastOfTie] / nctrl),
      y = c(0, cumsum(orderedLabels == 1)[lastOfTie] / ncase)
    )
    rocPlot <- ggplot(rocData, aes(x = x, y = y)) +
      geom_step(linewidth = 1) +
      geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 1.5) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_bw() +
      ylab("True Positive Rate (Sensitivity)") +
      xlab("False Positive Rate (1 - Specificity)")
  } else {
    auc     <- NULL
    rocPlot <- NULL
  }
  est       <- list(selection  = finalMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    proportion = proportion,
                    targetVec  = targetVector,
                    rocPlot    = rocPlot,
                    finalBest  = finalBestCrit,
                    optimize   = optimize)
  class(est) <- "calf"
  return(est)
}

compute.auc <- function(caseVar, ctrlVar){
  funcValue <- c(caseVar, ctrlVar)
  funcValue <- round(funcValue, digits = 8)
  ranks     <- rank(funcValue, ties.method = "average")
  caseFunc  <- sum(ranks[1:length(caseVar)]) - length(caseVar)*(length(caseVar)+1)/2
  ctrlFunc  <- sum(ranks[(length(caseVar)+1):length(ranks)]) - length(ctrlVar)*(length(ctrlVar)+1)/2
  auc       <- round(max(ctrlFunc, caseFunc)/(caseFunc + ctrlFunc), digits = 4)
  return(auc)
}

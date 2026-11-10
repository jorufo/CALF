validate_cv_inputs <- function(data, limit, proportion, times, targetVector, optimize, outputPath) {
  if (!is.data.frame(data) && !is.matrix(data))
    stop("CALF ERROR: data must be a matrix or data frame.")
  if (ncol(data) < 2L || nrow(data) < 6L ||
      any(!vapply(as.data.frame(data), is.numeric, logical(1))))
    stop("CALF ERROR: data must have numeric columns, at least one marker, and at least six rows.")
  if (length(limit) != 1L || !is.finite(limit) || limit != as.integer(limit) ||
      limit < 1L || limit > ncol(data) - 1L)
    stop("CALF ERROR: limit must be an integer between 1 and the number of markers.")
  if (length(times) != 1L || !is.finite(times) || times != as.integer(times) || times < 1L)
    stop("CALF ERROR: times must be a positive integer.")
  if (length(proportion) != 1L || !is.finite(proportion) || proportion <= 0 || proportion >= 1)
    stop("CALF ERROR: proportion must be strictly between zero and one.")
  if (length(targetVector) != 1L || !targetVector %in% c("binary", "nonbinary"))
    stop("CALF ERROR: targetVector must be 'binary' or 'nonbinary'.")
  allowed <- if (targetVector == "binary") c("pval", "auc") else "corr"
  if (length(optimize) != 1L || !optimize %in% allowed)
    stop(paste("CALF ERROR: optimize must be one of", paste(allowed, collapse = ", ")))
  if (targetVector == "binary") {
    if (anyNA(data[, 1]) || !all(data[, 1] %in% c(0, 1)))
      stop("CALF ERROR: binary target must contain only 0 and 1.")
    group_sizes <- c(sum(data[, 1] == 0), sum(data[, 1] == 1))
    kept <- round(group_sizes * proportion)
    if (any(kept < 2L) || any(group_sizes - kept < 2L))
      stop("CALF ERROR: proportion must leave at least two controls and cases in both sets.")
  } else {
    kept <- round(nrow(data) * proportion)
    if (kept < 3L || nrow(data) - kept < 3L)
      stop("CALF ERROR: proportion must leave at least three rows in both sets.")
  }
  if (!is.null(outputPath) &&
      (length(outputPath) != 1L || !dir.exists(outputPath)))
    stop("CALF ERROR: outputPath must be an existing directory.")
  invisible(TRUE)
}



#'@import data.table
#'@import ggplot2
#'@title calf
#'@description Coarse Approximation Linear Function
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "nonbinary"). All other columns contain relevant markers.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "nonbinary" for target vector with real numbers.
#'@param optimize Criteria to optimize, "pval" or "auc", (if targetVector = "binary") or "corr" (if targetVector = "nonbinary").  Defaults to "pval".
#'@param verbose Logical. Indicate TRUE to print activity at each iteration to console. Defaults to FALSE.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The optimal AUC, pval, or correlation for the classification.
#'@return If targetVector is binary, rocPlot. A plot object from ggplot2 for the receiver operating curve.
#'@examples
#'calf(data = CaseControl, nMarkers = 6, targetVector = "binary", optimize = "pval")
#'@export
calf <- function(data,
                 nMarkers,
                 targetVector,
                 optimize = "pval",
                 verbose = FALSE){
  calf_internal(data,
                nMarkers,
                proportion = NULL,
                randomize  = FALSE,
                targetVector = targetVector,
                times      = 1,
                optimize = optimize,
                verbose = verbose)
}



#'@title calf_fractional
#'@description Randomly selects from binary input provided to data parameter while ensuring the requested proportions of case and control variables are used and runs Coarse Approximation Linear Function.
#'@param data Matrix or data frame. Must be binary data such that the first column must contain case/control dummy coded variable, as function is only approprite for binary data.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param controlProportion Proportion of control samples to use, default is .8.
#'@param caseProportion Proportion of case samples to use, default is .8.
#'@param optimize Criteria to optimize, "pval" or "auc".  Defaults to "pval".
#'@param verbose Logical. Indicate TRUE to print activity at each iteration to console. Defaults to FALSE.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The optimal AUC or pval for the classification.
#'@return rocPlot. A plot object from ggplot2 for the receiver operating curve.
#'@examples
#'calf_fractional(data = CaseControl, nMarkers = 6, controlProportion = .8, caseProportion = .4)
#'@export
calf_fractional <- function(data,
                     nMarkers,
				             controlProportion = .8,
				             caseProportion = .8,
                     optimize = "pval",
                     verbose = FALSE){
  
     calf_internal(data,
                   nMarkers,
                   proportion = c(controlProportion,caseProportion),
                   randomize  = FALSE,
                   targetVector = "binary",
                   times      = 1,
                   optimize = optimize,
                   verbose = verbose)
}
				


#'@title calf_randomize
#'@description Randomly selects from binary input provided to data parameter and runs Coarse Approximation Linear Function.
#'@param data Matrix or data frame. Must be binary data such that the first column must contain case/control dummy coded variable, as function is only approprite for binary data.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "nonbinary" for target vector with real numbers.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param optimize Criteria to optimize if targetVector = "binary." Indicate "pval" to optimize the p-value corresponding to the t-test distinguishing case and control. Indicate "auc" to optimize the AUC.
#'@param verbose Logical. Indicate TRUE to print activity at each iteration to console. Defaults to FALSE.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The optimal AUC, pval, or correlation for the classification.
#'@return aucHist A histogram of the AUCs across replications, if applicable.
#'@examples
#'calf_randomize(data = CaseControl, nMarkers = 6, targetVector = "binary", times = 5)
#'@export
calf_randomize <- function(data,
                           nMarkers,
                           targetVector,
                           times      = 1,
                           optimize   = "pval",
                           verbose = FALSE){
  auc        <- numeric()
  finalBest  <- numeric()
  allMarkers <- character()
  count      <- 1
  AUC = NULL
  
  randomize = TRUE
  
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = NULL,
                         randomize  = randomize,
                         targetVector = targetVector,
                         times,
                         optimize = optimize,
                         verbose = verbose)
    if(!is.null(out$auc))
      auc[count] <- out$auc
    
    selection  <- out$selection
    markers    <- as.character(out$selection[,1])
    finalBest  <- append(finalBest, out$finalBest)
    allMarkers <- as.character((append(allMarkers, markers)))
    if (count == times) break
    count      <- count + 1
  }

  if (times > 1) {
    summaryMarkers <- as.data.frame(table(allMarkers), check.names = FALSE)
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    if (targetVector == "binary"){
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
    } else aucHist = NULL
  } else {
    summaryMarkers = NULL
    aucHist        = NULL
  }
  if (times == 1 & targetVector == "binary") {
    rocPlot <- out$rocPlot
  } else {
    rocPlot <- NULL
  }

  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    targetVec  = targetVector,
                    aucHist    = aucHist,
                    times      = times,
                    finalBest  = finalBest,
                    rocPlot    = rocPlot,
                    optimize   = optimize,
                    verbose    = verbose)
  class(est) <- "calf_randomize"
  return(est)
}



#'@title calf_subset
#'@description Runs Coarse Approximation Linear Function on a random subset of the data provided, resulting in the same proportion applied to case and control, when applicable.
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "nonbinary"). All other columns contain relevant markers.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param proportion Numeric. A value between 0 and 1 indicating the proportion of cases and controls to use in analysis (if targetVector = "binary"). If targetVector = "nonbinary", this is just a proportion of the full sample. Used to evaluate robustness of solution. Defaults to 0.8.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "nonbinary" for target vector with real numbers.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param optimize Criteria to optimize if targetVector = "binary." Indicate "pval" to optimize the p-value corresponding to the t-test distinguishing case and control. Indicate "auc" to optimize the AUC.
#'@param verbose Logical. Indicate TRUE to print activity at each iteration to console. Defaults to FALSE.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The optimal AUC, pval, or correlation for the classification. If multiple replications are requested, a data.frame containing all optimized values across all replications is returned.
#'@return aucHist A histogram of the AUCs across replications, if applicable.
#'@examples
#'calf_subset(data = CaseControl, nMarkers = 6, targetVector = "binary", times = 5)
#'@export

calf_subset <- function(data,
                        nMarkers,
                        proportion = .8,
                        targetVector,
                        times      = 1,
                        optimize = "pval",
                        verbose = FALSE){
  auc        <- numeric()
  allMarkers <- character()
  finalBest  <- numeric()
  count      <- 1
  AUC = NULL
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = proportion,
                         randomize  = FALSE,
                         targetVector = targetVector,
                         times,
                         optimize = optimize,
                         verbose = verbose)
    
    if(!is.null(out$auc))
      auc[count] <- out$auc

    selection  <- out$selection
    finalBest  <- append(finalBest, out$finalBest)
    markers    <- as.character(out$selection[,1])
    allMarkers <- as.character((append(allMarkers, markers)))
    if (count == times) break
    count      <- count + 1
  }

  if (times > 1){
    summaryMarkers <- as.data.frame(table(allMarkers), check.names = FALSE)
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    if (targetVector == "binary"){
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
    } else aucHist = NULL
  } else {
    summaryMarkers = NULL
    aucHist        = NULL
  }
  if (times == 1 & targetVector == "binary") {
    rocPlot <- out$rocPlot
  } else {
    rocPlot <- NULL
  }

  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    proportion = proportion,
                    targetVec  = targetVector,
                    aucHist    = aucHist,
                    times      = times,
                    finalBest  = finalBest,
                    rocPlot    = rocPlot,
                    optimize   = optimize)
  class(est) <- "calf_subset"
  return(est)
}







#'@title calf_exact_binary_subset
#'@description Runs Coarse Approximation Linear Function on a random subset of binary data provided, with the ability to precisely control the number of case and control data used.
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param nCase Numeric. A value indicating the number of case data to use.
#'@param nControl Numeric. A value indicating the number of control data to use.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param optimize Criteria to optimize.  Indicate "pval" to optimize the p-value corresponding to the t-test distinguishing case and control. Indicate "auc" to optimize the AUC.
#'@param verbose Logical. Indicate TRUE to print activity at each iteration to console. Defaults to FALSE.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The optimal AUC or pval for the classification. If multiple replications are requested, a data.frame containing all optimized values across all replications is returned.
#'@return aucHist A histogram of the AUCs across replications, if applicable.
#'@examples
#'calf_exact_binary_subset(data = CaseControl, nMarkers = 6, nCase = 5, nControl = 8, times = 5)
#'@export
calf_exact_binary_subset <- function(data,
                        nMarkers,
                        nCase,
                        nControl,
                        times      = 1,
                        optimize = "pval",
                        verbose = FALSE){
    if (length(nCase) != 1L || length(nControl) != 1L ||
        !is.finite(nCase) || !is.finite(nControl) ||
        nCase != as.integer(nCase) || nControl != as.integer(nControl) ||
        nCase < 2L || nControl < 2L)
      stop("CALF ERROR: nCase and nControl must be integers of at least two.")
    if (length(times) != 1L || !is.finite(times) ||
        times != as.integer(times) || times < 1L)
      stop("CALF ERROR: times must be a positive integer.")
  
    targetVector = "binary"
    proportion = 1
  
  
    #Determine which is case and which is control
    ctrlRows  <- which(data[ ,1] == 0)
    caseRows  <- which(data[ ,1] == 1)
    if (nControl > length(ctrlRows) || nCase > length(caseRows))
      stop("CALF ERROR: Requested cases or controls exceed available observations.")

    auc        <- numeric()
    allMarkers <- character()
    finalBest  <- numeric()
    count      <- 1
    AUC = NULL
    repeat {
      
      #Resample the binary data, thus controlling the randomization here.
      keepRows <- c(sample(ctrlRows, nControl), sample(caseRows, nCase))
      resampledData <- data[keepRows, , drop = FALSE]
      
      
      out <- calf_internal(resampledData,
                           nMarkers,
                           proportion = proportion,
                           randomize  = FALSE,
                           targetVector = targetVector,
                           times,
                           optimize = optimize,
                           verbose = verbose)
      auc[count] <- out$auc
      selection  <- out$selection
      finalBest  <- append(finalBest, out$finalBest)
      markers    <- as.character(out$selection[,1])
      allMarkers <- as.character((append(allMarkers, markers)))
      if (count == times) break
      count      <- count + 1
    }
    
    if (times > 1){
      summaryMarkers <- as.data.frame(table(allMarkers), check.names = FALSE)
      colnames(summaryMarkers) <- c("Marker", "Frequency")
      summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]

      auc            <- as.data.frame(auc)
      colnames(auc)  <- "AUC"
      aucHist <- ggplot(auc, aes(AUC)) +
        geom_histogram() +
        ylab("Count") +
        xlab("AUC") +
        scale_x_continuous() +
        theme_bw()

    } else {
      summaryMarkers = NULL
      aucHist        = NULL
    }
    
    if (times == 1) {
      rocPlot <- out$rocPlot
    } else {
      rocPlot <- NULL
    }
    
    est       <- list(selection  = selection,
                      multiple   = summaryMarkers,
                      auc        = auc,
                      proportion = proportion,
                      targetVec  = targetVector,
                      aucHist    = aucHist,
                      times      = times,
                      finalBest  = finalBest,
                      rocPlot    = rocPlot,
                      optimize   = optimize)
    class(est) <- "calf_exact_binary_subset"
    return(est)
}









#'@title cv.calf
#'@description Performs cross-validation using CALF data input
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "nonbinary"). All other columns contain relevant markers.
#'@param limit Maximum number of markers to include in creation of sum.
#'@param proportion Numeric. A value between 0 and 1 indicating the proportion of cases and controls to use in analysis (if targetVector = "binary") or proportion of the full sample (if targetVector = "nonbinary"). Defaults to 0.8.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "nonbinary" for target vector with real numbers.
#'@param optimize Criteria to optimize if targetVector = "binary." Indicate "pval" to optimize the p-value corresponding to the t-test distinguishing case and control. Indicate "auc" to optimize the AUC.  Defaults to pval.
#'@param outputPath The path where files are to be written as output, default is NULL meaning no files will be written.  When targetVector is "binary" file binary.csv will be output in the provided path, showing the reults.  When targetVector is "nonbinary" file nonbinary.csv will be output in the provided path, showing the results.  In the same path, the kept and unkept variables from the last iteration, will be output, prefixed with the targetVector type "binary" or "nonbinary" followed by Kept and Unkept and suffixed with .csv.  Two files containing the results from each run have List in the filenames and suffixed with .txt.
#'@return A data frame containing "times" rows of CALF runs where each row represents a run of CALF on a randomized "proportion" of "data".  Colunns start with the numer selected for the run, followed by AUC or pval and then all markers from "data".  An entry in a marker column signifys a chosen marker for a particular run (a row) and their assigned coarse weight (-1, 0, or 1).
#'@examples
#'\dontrun{
#'cv.calf(data = CaseControl, limit = 5, times = 100, targetVector = 'binary')
#'}
#'@export
cv.calf <- function(data, limit, proportion = .8, times, targetVector, optimize = "pval", outputPath=NULL) {
  validate_cv_inputs(data, limit, proportion, times, targetVector, optimize, outputPath)
  
  if (targetVector != "binary" && targetVector != "nonbinary") {
    cat('CALF ERROR: Invalid targetVector argument.  Only "binary" or "nonbinary" is allowed.')
  } else if (targetVector == "binary" && optimize=="corr") {
    cat('CALF ERROR: Optimizing by "corr" is only applicable to nonbinary data.')
  } else if (targetVector == "nonbinary" && optimize=="pval") {
    cat('CALF ERROR: Optimizing by "pval" is only applicable to binary data.')
  } else if (targetVector == "nonbinary" && optimize=="auc") {
    cat('CALF ERROR: Optimizing by "auc" is only applicable to binary data.')
  } else {
    
    #Get the rows of interest first as there is no reason to repeat this
    if (targetVector == "binary") {
      
      ctrlRows  <- which(data[ ,1] == 0)
      caseRows  <- which(data[ ,1] == 1)
      
      # calculate number of case and control to keep
      nCtrlKeep <- round(length(ctrlRows)*proportion, digits = 0)
      nCaseKeep <- round(length(caseRows)*proportion, digits = 0)
      
    } else if(targetVector == "nonbinary"){
      
      nDataKeep <- round(nrow(data)*proportion, digits = 0)
      
    } 
    

    #Build the header row for the table that will be output
    if (targetVector == "binary") {
      if (optimize == "pval") {
        header <- c("Number Selected", "AUC", "pval", colnames(data)[-1])
      } else if (optimize == "auc"){
        header <- c("Number Selected", "AUC", colnames(data)[-1])
      }
    } else if (targetVector == "nonbinary"){
      header <- c("Number Selected", "corr", colnames(data)[-1])
    }
    
    results <- matrix(0, times, length(header))
    colnames(results)<-header
    
    

    #Now run the CALF calculation "times" times
    rowCount = 1
    optimizedKeptList <- vector()
    optimizedUnkeptList <- vector()
    correlationList <- vector()
    repeat {
      
      if (targetVector == "binary") {
        
        #Resample the binary data, keeping track of what was included and what was not.
        keepCtrlRows <- sample(ctrlRows, nCtrlKeep)
        unkeptCtrlRows <- setdiff(union(ctrlRows,keepCtrlRows), intersect(ctrlRows,keepCtrlRows))
        
        keepCaseRows <- sample(caseRows, nCaseKeep)
        unkeptCaseRows <- setdiff(union(caseRows,keepCaseRows), intersect(caseRows,keepCaseRows))
        
        keepRows  <- c(keepCtrlRows, keepCaseRows)
        unkeptRows <- c(unkeptCtrlRows, unkeptCaseRows)
        
        unkeptCaseData <- data[unkeptCaseRows, , drop = FALSE]
        unkeptCtrlData <- data[unkeptCtrlRows, , drop = FALSE]
        
        resampledData <- data[keepRows, , drop = FALSE]
        unkeptData <- data[unkeptRows, , drop = FALSE]
        
        if(!is.null(outputPath)) {
          outputFile <- file.path(outputPath, "binaryKept.csv")
          fwrite(resampledData, outputFile)
          
          
          outputFile <- file.path(outputPath, "binaryUnkept.csv")
          fwrite(unkeptData, outputFile)
        }
        
      } else if(targetVector == "nonbinary"){
        
        #Resample the nonbinary data
        keepRows  <- sample(seq_len(nrow(data)), nDataKeep)
        unkeptRows <- setdiff(seq(1, length(data[,1]), by=1), keepRows)
        resampledData <- data[keepRows, , drop = FALSE]
        unkeptData <- data[unkeptRows, , drop = FALSE]
      
        if(!is.null(outputPath)) {
          outputFile <- file.path(outputPath, "nonbinaryKept.csv")
          fwrite(resampledData, outputFile)
          
          outputFile <- file.path(outputPath, "nonbinaryUnkept.csv")
          fwrite(unkeptData, outputFile)
        }
        
      }
      

      answer = calf_internal(data=resampledData,
                             nMarkers = limit,
                             randomize  = FALSE,
                             proportion = ,
                             times = 1,
                             targetVector = targetVector,
                             optimize = optimize,
                             verbose = FALSE)
      
      
      #Keep track of the optimizer values returned for each run
      if(optimize == "auc") {
        results[rowCount, "AUC"] = answer$auc
        optimizedKeptList <- append(optimizedKeptList, answer$auc)
      } else if(optimize == "pval") {
        results[rowCount, "AUC"] = answer$auc
        results[rowCount, "pval"] = answer$finalBest
        optimizedKeptList <- append(optimizedKeptList, answer$finalBest)
      } else if(optimize == "corr") {
        results[rowCount, "corr"] = answer$finalBest
        optimizedKeptList <- append(optimizedKeptList, answer$finalBest)
      }
      
      
      #Keep a tally of the results per calf run
      markerCount = 1
      markerList = as.character(answer$selection$Marker)
      lenMarkerList = length(markerList)
      results[rowCount, "Number Selected"] = lenMarkerList
      repeat {
        
        results[rowCount, markerList[markerCount]] = answer$selection$Weight[markerCount]
        
        markerCount <- markerCount + 1
        
        if (markerCount > lenMarkerList)
          break
      }
      
      

      #Perform the cross-validation
      if (targetVector == "binary") {
        if (optimize == "pval") {
          header <- c("Number Selected", "AUC", "pval", colnames(data)[-1])
          weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
            matrix(results[rowCount, -(1:3)], ncol = 1)
  
          resultCtrlData <- weightsTimesUnkept[seq_len(nrow(unkeptCtrlData))]
          resultCaseData <- weightsTimesUnkept[nrow(unkeptCtrlData) + seq_len(nrow(unkeptCaseData))]
  
          #optimizedUnkeptList<-append(optimizedUnkeptList, t.test(resultCaseData, resultCtrlData, var.equal = FALSE)$p.value)
          optimizedUnkeptList<-append(optimizedUnkeptList, compute.auc(resultCaseData, resultCtrlData))
          
        } else if (optimize == "auc"){
          weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
            matrix(results[rowCount, -(1:2)], ncol = 1)
      
          resultCtrlData <- weightsTimesUnkept[seq_len(nrow(unkeptCtrlData))]
          resultCaseData <- weightsTimesUnkept[nrow(unkeptCtrlData) + seq_len(nrow(unkeptCaseData))]
      
          optimizedUnkeptList<-append(optimizedUnkeptList, compute.auc(resultCaseData, resultCtrlData))
            
        }
      } else if (targetVector == "nonbinary"){
      
        weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
          matrix(results[rowCount, -(1:2)], ncol = 1)
        corrResult <- cor(as.vector(weightsTimesUnkept), unkeptData[[1]])
        correlationList <- append(correlationList,corrResult )
      }
    
    
      rowCount <- rowCount + 1
        
      if (rowCount > times)
        break
    }
  
  }
    
    
    
  #If an outputPath was provided, then output the extra data generated by the CV
  if(!is.null(outputPath)) {
    #Write the results
    if (targetVector == "binary") {
      
      outputFile <- file.path(outputPath, "binary.csv")
      fwrite(results, outputFile)
      
      outputFile <- file.path(outputPath, paste0(optimize, "KeptList.txt"))
      write(optimizedKeptList, outputFile, sep="\n")
      
      outputFile <- file.path(outputPath, "AUCUnkeptList.txt")
      write(optimizedUnkeptList, outputFile, sep="\n")
      
    } else if(targetVector == "nonbinary"){
      
      outputFile <- file.path(outputPath, "nonbinary.csv")
      fwrite(results, outputFile)
      
      outputFile <- file.path(outputPath, "corrUnkeptList.txt")
      write(correlationList, outputFile, sep="\n")
      
    }
    

    
  }


  return(results)
}







#'@title perm_target_cv.calf
#'@description Performs cross-validation using CALF data input and randomizes the target column with each iteration of the loop, controlled by 'times'.
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "nonbinary"). All other columns contain relevant markers.
#'@param limit Maximum number of markers to include in creation of sum.
#'@param proportion Numeric. A value between 0 and 1 indicating the proportion of cases and controls to use in analysis (if targetVector = "binary") or proportion of the full sample (if targetVector = "nonbinary"). Defaults to 0.8.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "nonbinary" for target vector with real numbers.
#'@param optimize Criteria to optimize if targetVector = "binary." Indicate "pval" to optimize the p-value corresponding to the t-test distinguishing case and control. Indicate "auc" to optimize the AUC.  Defaults to pval.
#'@param outputPath The path where files are to be written as output, default is NULL meaning no files will be written.  When targetVector is "binary" file binary.csv will be output in the provided path, showing the reults.  When targetVector is "nonbinary" file nonbinary.csv will be output in the provided path, showing the results.  In the same path, the kept and unkept variables from the last iteration, will be output, prefixed with the targetVector type "binary" or "nonbinary" followed by Kept and Unkept and suffixed with .csv.  Two files containing the results from each run have List in the filenames and suffixed with .txt.
#'@return A data frame containing "times" rows of CALF runs where each row represents a run of CALF on a randomized "proportion" of "data".  Colunns start with the numer selected for the run, followed by AUC or pval and then all markers from "data".  An entry in a marker column signifys a chosen marker for a particular run (a row) and their assigned coarse weight (-1, 0, or 1).
#'@examples
#'\dontrun{
#'perm_target_cv.calf(data = CaseControl, limit = 5, times = 100, targetVector = 'binary')
#'}
#'@export
perm_target_cv.calf <- function(data, limit, proportion = .8, times, targetVector, optimize = "pval", outputPath=NULL) {
  validate_cv_inputs(data, limit, proportion, times, targetVector, optimize, outputPath)
  
  if (targetVector != "binary" && targetVector != "nonbinary") {
    cat('CALF ERROR: Invalid targetVector argument.  Only "binary" or "nonbinary" is allowed.')
  } else if (targetVector == "binary" && optimize=="corr") {
    cat('CALF ERROR: Optimizing by "corr" is only applicable to nonbinary data.')
  } else if (targetVector == "nonbinary" && optimize=="pval") {
    cat('CALF ERROR: Optimizing by "pval" is only applicable to binary data.')
  } else if (targetVector == "nonbinary" && optimize=="auc") {
    cat('CALF ERROR: Optimizing by "auc" is only applicable to binary data.')
  } else {
    
    #Get the rows of interest first as there is no reason to repeat this
    if (targetVector == "binary") {
      
      ctrlRows  <- which(data[ ,1] == 0)
      caseRows  <- which(data[ ,1] == 1)
      
      # calculate number of case and control to keep
      nCtrlKeep <- round(length(ctrlRows)*proportion, digits = 0)
      nCaseKeep <- round(length(caseRows)*proportion, digits = 0)
      
    } else if(targetVector == "nonbinary"){
      
      nDataKeep <- round(nrow(data)*proportion, digits = 0)
      
    } 
    
    
    #Build the header row for the table that will be output
    if (targetVector == "binary") {
      if (optimize == "pval") {
        header <- c("Number Selected", "AUC", "pval", colnames(data)[-1])
      } else if (optimize == "auc"){
        header <- c("Number Selected", "AUC", colnames(data)[-1])
      }
    } else if (targetVector == "nonbinary"){
      header <- c("Number Selected", "corr", colnames(data)[-1])
    }
    
    results <- matrix(0, times, length(header))
    colnames(results)<-header
    
    
    
    #Now run the CALF calculation "times" times
    rowCount = 1
    optimizedKeptList <- vector()
    optimizedUnkeptList <- vector()
    correlationList <- vector()
    repeat {
      
      print(paste("Iteration: ", rowCount))
      
      if (targetVector == "binary") {
        
        #Resample the binary data, keeping track of what was included and what was not.
        
        keepCtrlRows <- sample(ctrlRows, nCtrlKeep)
        unkeptCtrlRows <- setdiff(ctrlRows, keepCtrlRows)
        
        
        keepCaseRows <- sample(caseRows, nCaseKeep)
        unkeptCaseRows <- setdiff(caseRows, keepCaseRows)
        
        keepRows  <- c(keepCtrlRows, keepCaseRows)
        unkeptRows <- c(unkeptCtrlRows, unkeptCaseRows)
        
        unkeptCaseData <- data[unkeptCaseRows, , drop = FALSE]
        unkeptCtrlData <- data[unkeptCtrlRows, , drop = FALSE]
        
        resampledData <- data[keepRows, , drop = FALSE]
        unkeptData <- data[unkeptRows, , drop = FALSE]
        # Permute training labels; holdout labels remain the ground truth.
        resampledData[, 1] <- sample(resampledData[, 1])
        
        if(!is.null(outputPath)) {
          outputFile <- file.path(outputPath, "binaryKept.csv")
          fwrite(resampledData, outputFile)
          
          
          outputFile <- file.path(outputPath, "binaryUnkept.csv")
          fwrite(unkeptData, outputFile)
        }
        
      } else if(targetVector == "nonbinary"){
        
        keepRows  <- sample(seq_len(nrow(data)), nDataKeep)
        unkeptRows <- setdiff(seq_len(nrow(data)), keepRows)

        resampledData <- data[keepRows, , drop = FALSE]
        unkeptData <- data[unkeptRows, , drop = FALSE]
        resampledData[, 1] <- sample(resampledData[, 1])
        
        if(!is.null(outputPath)) {
          outputFile <- file.path(outputPath, "nonbinaryKept.csv")
          fwrite(resampledData, outputFile)
          
          outputFile <- file.path(outputPath, "nonbinaryUnkept.csv")
          fwrite(unkeptData, outputFile)
        }
        
      }
      
      
      answer = calf_internal(data=resampledData,
                             nMarkers = limit,
                             randomize  = FALSE,
                             proportion = ,
                             times = 1,
                             targetVector = targetVector,
                             optimize = optimize,
                             verbose = FALSE)
      
      
      #Keep track of the optimizer values returned for each run
      if(optimize == "auc") {
        results[rowCount, "AUC"] = answer$auc
        optimizedKeptList <- append(optimizedKeptList, answer$auc)
      } else if(optimize == "pval") {
        results[rowCount, "AUC"] = answer$auc
        results[rowCount, "pval"] = answer$finalBest
        optimizedKeptList <- append(optimizedKeptList, answer$finalBest)
      } else if(optimize == "corr") {
        results[rowCount, "corr"] = answer$finalBest
        optimizedKeptList <- append(optimizedKeptList, answer$finalBest)
      }
      
      
      #Keep a tally of the results per calf run
      markerCount = 1
      markerList = as.character(answer$selection$Marker)
      lenMarkerList = length(markerList)
      results[rowCount, "Number Selected"] = lenMarkerList
      repeat {
        
        results[rowCount, markerList[markerCount]] = answer$selection$Weight[markerCount]
        
        markerCount <- markerCount + 1
        
        if (markerCount > lenMarkerList)
          break
      }
      
      
      
      #Perform the cross-validation
      if (targetVector == "binary") {
        if (optimize == "pval") {
          header <- c("Number Selected", "AUC", "pval", colnames(data)[-1])
          weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
            matrix(results[rowCount, -(1:3)], ncol = 1)
          
          resultCtrlData <- weightsTimesUnkept[seq_len(nrow(unkeptCtrlData))]
          resultCaseData <- weightsTimesUnkept[nrow(unkeptCtrlData) + seq_len(nrow(unkeptCaseData))]
          
          #optimizedUnkeptList<-append(optimizedUnkeptList, t.test(resultCaseData, resultCtrlData, var.equal = FALSE)$p.value)
          optimizedUnkeptList<-append(optimizedUnkeptList, compute.auc(resultCaseData, resultCtrlData))
          
        } else if (optimize == "auc"){
          weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
            matrix(results[rowCount, -(1:2)], ncol = 1)
          
          resultCtrlData <- weightsTimesUnkept[seq_len(nrow(unkeptCtrlData))]
          resultCaseData <- weightsTimesUnkept[nrow(unkeptCtrlData) + seq_len(nrow(unkeptCaseData))]
          
          optimizedUnkeptList<-append(optimizedUnkeptList, compute.auc(resultCaseData, resultCtrlData))
          
        }
      } else if (targetVector == "nonbinary"){
        
        weightsTimesUnkept <- as.matrix(unkeptData[, -1, drop = FALSE]) %*%
          matrix(results[rowCount, -(1:2)], ncol = 1)
        corrResult <- cor(as.vector(weightsTimesUnkept), unkeptData[[1]])
        correlationList <- append(correlationList,corrResult )
      }
      
      
      rowCount <- rowCount + 1
      
      if (rowCount > times)
        break
    }
    
  }
  
  
  
  #If an outputPath was provided, then output the extra data generated by the CV
  if(!is.null(outputPath)) {
    #Write the results
    if (targetVector == "binary") {
      
      outputFile <- file.path(outputPath, "binary.csv")
      fwrite(results, outputFile)
      
      outputFile <- file.path(outputPath, paste0(optimize, "KeptList.txt"))
      write(optimizedKeptList, outputFile, sep="\n")
      
      outputFile <- file.path(outputPath, "AUCUnkeptList.txt")
      write(optimizedUnkeptList, outputFile, sep="\n")
      
    } else if(targetVector == "nonbinary"){
      
      outputFile <- file.path(outputPath, "nonbinary.csv")
      fwrite(results, outputFile)
      
      outputFile <- file.path(outputPath, "corrUnkeptList.txt")
      write(correlationList, outputFile, sep="\n")
      
    }
    
    
    
  }
  
  
  return(results)
}




#'@title write.calf
#'@description Writes output of the CALF dataframe 
#'@param x A CALF data frame.
#'@param filename The output filename
#'@export
write.calf <- function(x, filename){
  
  write.table(x$selection,
              file = filename,
              sep = ",",
              row.names = FALSE)
  
  
  if(x$targetVec == "binary" && x$optimize=="auc") {
    
    write( paste("\n","AUC ,",x$finalBest),
           file = filename,
           append = TRUE)
    
  } else if(x$targetVec == "binary" && x$optimize=="pval") {
    
    write( paste("\n","pval ,",x$finalBest),
           file = filename,
           append = TRUE)
    
  } else if(x$targetVec == "nonbinary") {
    
    write( paste("\n","corr,", x$finalBest),
           file = filename,
           append = TRUE)
  }
  
  
  
  
}





#'@title write.calf_randomize
#'@description Writes output of the CALF randomize dataframe 
#'@param x A CALF randomize data frame.
#'@param filename The output filename
#'@export
write.calf_randomize <- function(x, filename){
  oldWarn <- getOption("warn")
  on.exit(options(warn = oldWarn), add = TRUE)
  options(warn=-1)
  
  write.table(x$selection,
              file = filename,
              sep = ",",
              row.names = FALSE)
  
  write("\n",
        file = filename,
        append = TRUE)
  
  write.table(x$multiple,
              file = filename,
              sep = ",",
              row.names = FALSE,
              append = TRUE)
  
  write("\n",
        file = filename,
        append = TRUE)
  
  if(x$targetVec == "binary" && x$optimize=="auc") {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("AUC")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
    
  } else if(x$targetVec == "binary" && x$optimize=="pval") {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("pval")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
    
  } else if(x$targetVec == "nonbinary") {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("corr")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
  }
  
}





#'@title write.calf_subset
#'@description Writes output of the CALF subset dataframe 
#'@param x A CALF subset data frame.
#'@param filename The output filename
#'@export
write.calf_subset <- function(x, filename){
  oldWarn <- getOption("warn")
  on.exit(options(warn = oldWarn), add = TRUE)
  options(warn=-1)
  
  write.table(x$selection,
              file = filename,
              sep = ",",
              row.names = FALSE)
  
  write("\n",
        file = filename,
        append = TRUE)
  
  write.table(x$multiple,
              file = filename,
              sep = ",",
              row.names = FALSE,
              append = TRUE)
  
  write("\n",
        file = filename,
        append = TRUE)
  
  if(x$targetVec == "binary" && (x$optimize=="auc")) {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("AUC")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
    
  } else if(x$targetVec == "binary" && x$optimize=="pval") {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("pval")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
    
  } else if(x$targetVec == "nonbinary") {
    
    finalBest = as.data.frame(x$finalBest)
    
    colnames(finalBest) <- c("corr")
    
    write.table( finalBest,
                 file = filename,
                 sep = ",",
                 append = TRUE)
  }
  
}


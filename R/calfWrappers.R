#' @title calf
#' @description Coarse Approximation Linear Function
#' 
#' @param resp Response vector; binary or continuous.
#' @param pred Matrix or data.frame with predictor variables.
#' @param maxPred Maximum number of predictors to include in the model.
#' @param score Either "auc", "pcc", or "wts"; see details for more information.
#' @param margin Real number from 0 to 1, the required improvement in score 
#' include an additional predictor; see details.
#' @param verbose Logical. When TRUE activity is printed for each iteration.
#'
#' @return 
#' A calf object with the following list items:
#' \describe{
#'   \item{\code{predInd}}{Index of selected predictors in \code{pred} in the 
#'   order of selection}
#'   \item{\code{predName}}{Name of the selected predictors in \code{pred} in
#'   the order of selection}
#'   \item{\code{weight}}{The weight value for the selected predictors in 
#'   the order of selection}
#'   \item{\code{trainVec}}{The best score value at each iteration with the
#'   addition of the new predictor given in predInd; included to monitor
#'   the training progress}
#'   \item{\code{score}}{The score function used to calculate the scores in
#'   \code{trainVec}}
#'   \item{\code{maxPred}}{The maximum number of predictors allowed to be 
#'   included in the model}
#'   \item{\code{margin}}{The required improvement in score to continue 
#'   adding parameters to the model}
#' }
#'
#' @details
#' 
#' \code{resp} must be logical for binary score functions, i.e. 
#' \code{score = "auc"} or \code{score = "wts"}. Otherwise, \code{resp} is 
#' considered continuous. 
#' 
#' Need to add lots here... 
#'
#' @examples
#' 
#' data(CalfBin)
#' calf(resp = as.logical(CalfBin[[1]]), pred = CalfBin[ , -1], maxPred = 5)
#'
#' @export

calf <- function(resp,
                 pred,
                 maxPred = ncol(pred),
                 score = "pcc",
                 margin = NULL,
                 verbose = FALSE) {
  
  calculateCalf(resp = resp,
                pred = pred,
                maxPred = maxPred,
                score = score,
                margin = margin,
                verbose = verbose)
  
}



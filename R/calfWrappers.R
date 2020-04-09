#' @title calf
#' @description Coarse Approximation Linear Function
#' 
#' @param x Matrix or data.frame with predictor variables.
#' @param y Response vector; binary or continuous.
#' @param maxX Maximum number of predictors to include in the model.
#' @param score Either "auc", "pcc", or "wts"; see details for more information.
#' @param margin Real number from 0 to 1, the required improvement in score 
#' include an additional predictor; see details.
#' @param verbose Logical. When TRUE activity is printed for each iteration.
#'
#' @return 
#' A calf object with the following list items:
#' \describe{
#'   \item{\code{xInd}}{Index of selected predictors in \code{x} in the 
#'   order of selection}
#'   \item{\code{xName}}{Name of the selected predictors in \code{x} in
#'   the order of selection}
#'   \item{\code{weight}}{The weight value for the selected predictors in 
#'   the order of selection}
#'   \item{\code{finalS}}{The final best score value at the end of optimization}
#'   \item{\code{trainVec}}{The best score value at each iteration with the
#'   addition of the new predictor given in predInd; included to monitor
#'   the training progress}
#'   \item{\code{score}}{The score function used to calculate the scores in
#'   \code{trainVec}}
#'   \item{\code{maxX}}{The maximum number of predictors allowed to be 
#'   included in the model}
#'   \item{\code{margin}}{The required improvement in score to continue 
#'   adding parameters to the model}
#' }
#'
#' @details
#' 
#' \code{y} must be logical for binary score functions, i.e. 
#' \code{score = "auc"} or \code{score = "wts"}. Otherwise, \code{y} is 
#' considered continuous. 
#' 
#' Need to add lots here... 
#'
#' @examples
#' 
#' data(CalfBin)
#' calf(x = CalfBin[ , -1], y = as.logical(CalfBin[[1]]), maxX = 5)
#'
#' @export

calf <- function(x,
                 y,
                 maxX = ncol(x),
                 score = "pcc",
                 margin = NULL,
                 verbose = FALSE) {
  
  calculateCalf(x = x,
                y = y,
                maxX = maxX,
                score = score,
                margin = margin,
                verbose = verbose)
  
}



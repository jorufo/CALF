#' @name scoreFuncs
#' @title CALF score functions
#' @description Available CALF score functions
#' @param resp Response vector
#' @param pred Predictor vector
#' @importFrom stats t.test cor
#' 
#' @details
#' \code{resp} MUST be logical for \code{calfScoreAUC} and \code{calfScoreWTS}.
#' 
#' Can add more later...
#' 

NULL

#' @describeIn scoreFuncs Calculate score as the area under the receivor 
#' operating curve; only for binary response
#' @export

calfScoreAUC <- function(resp, pred) {
  ranks <- rank(pred, ties.method = "average")
  n1 <- sum(resp)
  n2 <- sum(!resp)
  s1 <- sum(ranks[resp]) - n1*(n1 + 1)/2
  s2 <- sum(ranks[!resp]) - n2*(n2 + 1)/2
  s1/(s1 + s2)
}

#' @describeIn scoreFuncs Calculate score as the p-value from Welch's t-test; 
#' only for binary response
#' @export

calfScoreWTS <- function(resp, pred)  {
  -1*log(t.test(pred ~ resp)$p.value)
}

#' @describeIn scoreFuncs Calculate score as the absolute value of Pearson's 
#' correlation coefficient
#' @export

calfScorePCC <- function(resp, pred)  {
  abs(cor(pred, resp, use = "complete.obs"))
}
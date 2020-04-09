#' @name scoreFuncs
#' @title CALF score functions
#' @description Available CALF score functions
#' @param x Predictor vector
#' @param y Response vector
#' @importFrom stats t.test cor
#' 
#' @details
#' \code{y} MUST be logical for \code{calfScoreAUC} and \code{calfScoreWTS}.
#' 
#' Can add more later...
#' 

NULL

#' @describeIn scoreFuncs Calculate score as the area under the receivor 
#' operating curve; only for binary response
#' @export

calfScoreAUC <- function(x, y) {
  ranks <- rank(x, ties.method = "average")
  n1 <- sum(y)
  n2 <- sum(!y)
  s1 <- sum(ranks[y]) - n1*(n1 + 1)/2
  s2 <- sum(ranks[!y]) - n2*(n2 + 1)/2
  s1/(s1 + s2)
}

#' @describeIn scoreFuncs Calculate score as the p-value from Welch's t-test; 
#' only for binary response
#' @export

calfScoreWTS <- function(x, y)  {
  -1*log(t.test(x ~ y)$p.value)
}

#' @describeIn scoreFuncs Calculate score as the absolute value of Pearson's 
#' correlation coefficient
#' @export

calfScorePCC <- function(x, y)  {
  abs(cor(x, y, use = "complete.obs"))
}
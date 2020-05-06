#' @name CALF-package
#' @aliases CALF-package
#' @title Coarse Approximation Linear Function
#' @description Forward selection linear regression greedy algorithm.
#' @author {John Ford [aut] \cr,
#'    Stephanie Lane [aut, cre],\cr
#'    Clark Jeffries [aut], \cr
#'    Diana Perkins [aut]
#' }
#' Maintainer: John Ford \email{JoRuFo@@gmail.com}
#' @importFrom stats t.test cor
#' @importFrom utils write.table
#' @import ggplot2
#' @keywords calf
#' @details The Coarse Approximation Linear Function (CALF) algorithm is a type of forward selection
#' linear regression greedy algorithm. Nonzero weights are restricted to the values +1 and -1 and
#' their number limited by a parameter. Samples are controls (at least 2) and cases (at least 2).
#' A data matrix consists of a distinguished column that labels every row as either a control (0) or a case (1).
#' Other columns (at least one) contain real number marker measurement data.
#' Another input is a limit (positive integer) on the number of markers that can be selected for use in a linear sum.
#' The present version uses as a score of differentiation the two-tailed, two sample unequal variance Student t-test p-value.
#' Thus, any real-valued function applied to all samples generates values for controls and cases that are used to calculate the score.
#' CALF selects the one marker (first in case of tie) that best distinguishes controls from cases (score is smallest p-value).
#' CALF then checks the limit. If the number of selected markers is the limit, CALF ends.
#' Else, CALF seeks a second marker, if any, that best improves the score of the sum function generated
#' by adding the newly selected marker to the previous markers with weight +1 or weight -1.
#' The process continues until the limit is reached or until no additional marker can be included in the sum to improve the score.
#' The proportion parameter can take either a single value or a two column vector of values.  For a single value, the proportion supplied
#' is applied equally to the set of control data and the set of case data, and that desired propotion of data will be randomly selected 
#' from both the control data set and the case data set.  If a two column vector is supplied the first column value will be used as the 
#' proportion of the control data to use where the second value will represent the proportion of the case data to use.  The requested
#' propotions of data will be randomly selected from both the control data set and the case data set.
NULL

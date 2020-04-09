#' @docType data
#' @keywords datasets
#' @name CalfBin
#' @title Example binary data for CALF algorithm
#'
#' @description 
#' This data contains standardized data for 135 blood markers for 72 individuals 
#' who are distinguished as control (\code{CalfBin$NCvC == 0}) & case 
#' (\code{CalfBin$NCvC == 0}).
#' 
#' The response variable (designating case/control) is given by the first 
#' column, \code{NCvC}.The remaining columns, (i.e. \code{M---}) are different
#' blood markers (the predictor variables).
#' 
#' This data comes from the 
#' \href{https://napls.ucsf.edu}{North American Prodrome Longitudinal Study}.
#'
#' @usage data(CalfBin)
#' @format A data.table with the first column giving the response variable, 
#' and the remaining columns giving predictor variables.

NULL

#' @docType data
#' @keywords datasets
#' @name CalfCnt
#' @title Example continuous data for CALF algorithm
#'
#' @description  
#' This data contains standardized expression for 132 miRNAs in blood for 24 
#' people and their percentage of blood monocytes.
#' 
#' The response variable (fraction of monocytes) is given by the first 
#' column, \code{mono}.The remaining columns, (i.e. \code{M---}) are different
#' blood miRNAs (the predictor variables).
#' 
#' This data comes from the 
#' \href{https://napls.ucsf.edu}{North American Prodrome Longitudinal Study}.
#'
#' @usage data(CalfCnt)
#' @format A data.table with the first column giving the response variable, 
#' and the remaining columns giving predictor variables.

NULL

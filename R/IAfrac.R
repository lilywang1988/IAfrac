#' IAfrac: A package for interim analysis using weighted log-rank tests
#'
#' The IAfrac package implements Hasegawa (2014) and (2016) proposals for weighted log-rank tests in piece-wise expornential distributed survival functions. The current version only considers two pieces with only one change point (epsilon). This R package provides four categories of important functions:
#' sample size, information fraction, weighted log-rank test, data manipulation.
#' @section Sample size functions:
#' The sample size calculation functions implement methods proposed in Hasegawa (2014). They are recorded in R/Hasegama2014.R.
#'
#' @section Information fraction functions:
#' The iformaiton fraction functions implement methods proposed in Hasegawa (2016). They are recorded in R/Hasegama2016.R.
#'
#' @section  Weightd log-rank test:
#' The weighted log-rank test functions are newly written functions for weighted log-rank tests. They are recorded in R/WLRT.R.
#'
#' @section Data manipulation functions:
#' The data manipulation functions are prepared to trim the data either according to the follow-up time or the event counts. They are recorded in R/data_manipulation.R.
#' @docType package
#' @name IAfrac
NULL
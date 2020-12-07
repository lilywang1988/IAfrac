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
#' @import data.table
#' @import gsDesign
#' @export data.trim 
#' @export data.trim.d
#' @export sample.size_FH
#' @export avg.haz
#' @export getc
#' @export h.tilde
#' @export uv
#' @export v
#' @export u
#' @export h1
#' @export h0
#' @export I.0
#' @export I.0.cov 
#' @export I.1 
#' @export I.1.cov
#' @export cor.0
#' @export I_t
#' @export I_t.2
#' @export FH.frac.cal
#' @export WLR.test.cov
#' @export WLR.test.cor 
#' @export approx.I
#' @export FH.table
#' @export logrank.table
#' @export survKM_minus
#' @export survKM_exact
#' @export FH.test
#' @export WLR.test
#' 
NULL

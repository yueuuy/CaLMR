#' Example correlation matrix input of the function: samplecorr
#'
#' A 7*7 correlation matrix for six traits (B1 to B6) and exposure (Y).
#' The first six rows and columns represent the pairwise correlations between the traits.
#' The last row and column represent the correlations between each trait and the exposure (Y).
#'
#' @format A numeric matrix with 7 rows and 7 columns:
#' \describe{
#'   \item{B1 to B6}{Observable Traits (columns 1-6 and rows 1-6)}
#'   \item{Y}{The exposure variable (column 7 and row 7).}
#' }
#' @source Simulated correlation matrix for the package example.
"samplecorr"
